"""Exercise the REAL driver/core/panel with only external tool execution mocked.

Synthetic outputs test orchestration and checksums, NOT biological accuracy.
"""
import configparser
import contextlib
import csv
import importlib.util
import io
import json
import os
from pathlib import Path
import signal
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'pipeline'))
import pipeline as core
import filter_engine as engine
spec = importlib.util.spec_from_file_location('mock_gene_driver', ROOT/'run_gene_annotation.py')
driver = importlib.util.module_from_spec(spec); spec.loader.exec_module(driver)


def summary(path, lineage, c):
    path.mkdir(parents=True, exist_ok=True)
    (path/'summary.txt').write_text(f'lineage: {lineage}_odb12.2\nS: {c:.2f}%, {c}\nD: 0.00%, 0\nF: 0.00%, 0\nM: {100-c:.2f}%, {100-c}\nI: 0.00%, 0\nN: 100\n')
    (path/'full_table.tsv').write_text(''.join(f'm{i}\t'+('Single' if i<c else 'Missing')+'\n' for i in range(100)))


def annotation(dest, genes):
    dest.mkdir(parents=True, exist_ok=True)
    gtf, aa, cds, support = [], [], [], ['transcript_id\tnum_introns\tintrons_sup_any\n']
    for i in range(genes):
        attrs = f'gene_id "g{i}"; transcript_id "g{i}.t1";'
        for feature in ['gene', 'transcript', 'CDS']:
            gtf.append(f'chr1\tsynthetic\t{feature}\t{i*100+1}\t{i*100+90}\t.\t+\t0\t{attrs}\n')
        aa.append(f'>g{i}.t1\nMAAA\n'); cds.append(f'>g{i}.t1\nATGGCTGCTGCT\n')
        support.append(f'g{i}.t1\t1\t1\n')
    for name, text in [('braker.gtf', ''.join(gtf)), ('braker.aa', ''.join(aa)),
                       ('braker.codingseq', ''.join(cds)), ('braker.gff3', '##gff-version 3\n'),
                       ('gene_support.tsv', ''.join(support))]: (dest/name).write_text(text)


@unittest.skipUnless(sys.platform.startswith('linux'), 'Requires POSIX affinity/locks/safe removal')
class MockWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='gene_mock_'); self.root = Path(self.tmp.name).resolve()
        self.affinity = os.sched_getaffinity(0)
        self.handlers = {s:signal.getsignal(s) for s in (signal.SIGINT, signal.SIGTERM)}
        self.calls, self.fastp = [], []
        self.exe = self.file('bin/tool', '#!/bin/sh\nexit 0\n'); self.exe.chmod(0o755)
        self.images = {k:str(self.file('images/'+k+'.sif')) for k in ['braker3_image','agat_image','busco_image']}
        db = self.root/'db'
        for lin in ['arachnida','arthropoda']:
            for name in ['lengths_cutoff','scores_cutoff','refseq_db.faa.gz','hmms/test.hmm']:
                self.file('db/'+lin+'_odb12.2/'+name)
        self.cfg = dict(output_root=str(self.root/'output'), work_root=str(self.root/'scratch'), threads=1,
            memory_mb=64, rule_memory_mb=32, fastp_threads=1, cleanup_work=True, max_consecutive_failures=3,
            min_free_work_gb=0, min_free_output_gb=0, tools={k:str(self.exe) for k in ['snakemake','apptainer','compleasm','fastp']},
            containers=self.images, databases=dict(compleasm=str(db), busco=str(db)),
            container_cache=str(self.root/'cache/snakemake'), apptainer_cache=str(self.root/'cache/apptainer'), java_home='')
    def tearDown(self):
        os.sched_setaffinity(0,self.affinity)
        for signum, handler in self.handlers.items(): signal.signal(signum,handler)
        self.tmp.cleanup()
    def file(self, name, text='mock\n'):
        p=self.root/name; p.parent.mkdir(parents=True,exist_ok=True); p.write_text(text); return p
    def sample(self, name, mode='ep'):
        g=self.file('input/'+name+'.fa','>chr1\nACGTacgt\n'); p=self.file('input/'+name+'.aa','>p1\nMAAA\n')
        target=self.root/'input'/name/'arachnida'; summary(target,'arachnida',95)
        sample=dict(sample_id=name,mode=mode,genome=str(g),protein_fasta=str(p),genome_compleasm_summary=str(target/'summary.txt'),rnaseq_manifest='-',runs=[])
        if mode=='etp':
            runs=[]
            for n in range(2): runs.append(dict(run_id=f'run{n}',fastq_r1=str(self.file(f'input/{name}_{n}_1.fq')),fastq_r2=str(self.file(f'input/{name}_{n}_2.fq'))))
            sample.update(runs=runs,rnaseq_manifest=str(self.file('input/'+name+'.rna.tsv')))
        return sample
    def external(self, cmd, log, env=None, cwd=None):
        cmd=list(map(str,cmd)); self.calls.append(cmd)
        Path(log).parent.mkdir(parents=True,exist_ok=True); Path(log).write_text('mock external command\n')
        if '--dry-run' in cmd: return
        if 'protein' in cmd:
            lin=cmd[cmd.index('-l')+1]; out=Path(cmd[cmd.index('-o')+1])
            label=out.parent.name
            c=96 if label=='F' else 95 if label=='D' else 90
            summary(out,lin,c); return
        if '--snakefile' in cmd:
            work=Path(cmd[cmd.index('--directory')+1]); ini=configparser.ConfigParser(); ini.read(env['BRAKER4_CONFIG'])
            with open(ini['paths']['samples_file']) as stream: sample=next(csv.DictReader(stream))['sample_name']
            dest=work/'output'/sample
            annotation(dest, 3 if work.name=='F' else 5)
            if work.name=='01_braker':
                for rel in engine.INPUTS:
                    p=dest/rel
                    if not p.exists(): p.parent.mkdir(parents=True,exist_ok=True); p.write_text('mock upstream\n')
                (dest/'tsebra.log').write_text('single-exon gene filter active\n')
                (dest/'results').mkdir(); (dest/'results/.done').touch()
            return
        if '--out' in cmd:
            out=Path(cmd[cmd.index('--out')+1]); out.write_text('gene_id "g1"; transcript_id "g1.t1";\n'); return
        self.fail('Unexpected command '+repr(cmd))
    def fake_fastp(self, argv, log, environment):
        self.fastp.append(list(map(str,argv)))
        for flag in ['--out1','--out2','--html']: Path(argv[argv.index(flag)+1]).write_text('mock\n')
        Path(argv[argv.index('--json')+1]).write_text(json.dumps(dict(summary=dict(after_filtering=dict(total_reads=10)))))
    def execute(self, selected):
        with patch.object(engine,'call',side_effect=self.external), patch.object(driver,'native_call',side_effect=self.fake_fastp), contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return driver.run(self.cfg,selected,{})
    def test_ep_full_panel_skip_and_tamper(self):
        s=self.sample('sample_ep'); self.assertEqual(self.execute([s]),0)
        out=Path(self.cfg['output_root'])/'results/sample_ep'
        self.assertEqual(core.verify_completed(out)['selected'],'F')
        self.assertFalse(Path(self.cfg['work_root']).exists()); self.assertEqual(self.fastp,[])
        n=len(self.calls); self.assertEqual(self.execute([s]),0); self.assertEqual(len(self.calls),n)
        (out/'04_final/selected/braker.aa').write_text('tampered\n')
        self.assertEqual(self.execute([s]),1)
    def test_etp_pairs_and_serial_mixed_modes(self):
        selected=[self.sample('sample_etp','etp'),self.sample('sample_ep','ep')]
        self.assertEqual(self.execute(selected),0); self.assertEqual(len(self.fastp),2)
        for cmd in self.fastp:
            for arg in driver.guards.FASTP_ARGS: self.assertIn(arg,cmd)
        for s in selected: core.verify_completed(Path(self.cfg['output_root'])/'results'/s['sample_id'])
    def test_three_failures_stop_batch(self):
        selected=[self.sample('s'+str(i)) for i in range(4)]
        with patch.object(driver,'run_sample',side_effect=RuntimeError('synthetic failure')) as mocked, contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(driver.run(self.cfg,selected,{}),3); self.assertEqual(mocked.call_count,3)
    def test_single_failure_continues(self):
        selected=[self.sample('s'+str(i)) for i in range(3)]
        with patch.object(driver,'run_sample',side_effect=[RuntimeError('fail'),'completed','completed']) as mocked, contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(driver.run(self.cfg,selected,{}),1); self.assertEqual(mocked.call_count,3)
    def test_skipped_does_not_reset_failure_streak(self):
        selected=[self.sample('s'+str(i)) for i in range(5)]
        with patch.object(driver,'run_sample',side_effect=[RuntimeError('fail'),'skipped',RuntimeError('fail'),RuntimeError('fail'),'completed']) as mocked, contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(driver.run(self.cfg,selected,{}),3); self.assertEqual(mocked.call_count,4)


if __name__=='__main__': unittest.main()
