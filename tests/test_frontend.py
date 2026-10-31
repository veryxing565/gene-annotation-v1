"""Manifest validation, case preservation, ownership and resource tests."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('gene_driver', ROOT / 'run_gene_annotation.py')
driver = importlib.util.module_from_spec(spec); spec.loader.exec_module(driver)
guards = driver.guards


@unittest.skipUnless(sys.platform.startswith('linux'), 'Supported runtime is Linux')
class FrontendTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='gene_frontend_')
        self.root = Path(self.tmp.name).resolve()
    def tearDown(self): self.tmp.cleanup()
    def file(self, name, text='x\n'):
        p = self.root / name; p.parent.mkdir(parents=True, exist_ok=True); p.write_text(text); return p
    def manifest(self, mode='ep', rnaseq='-'):
        genome = self.file('g.fa', '>chr1\nACGTTacgtN\n')
        protein = self.file('p.fa', '>p1\nMAAA\n')
        summary = self.file('summary.txt', 'lineage: arachnida_odb12.2\n')
        return self.file('samples.tsv', '\t'.join(guards.SAMPLE_COLUMNS)+'\n'+
                         '\t'.join(['sample1', mode, str(genome), str(protein), str(summary), str(rnaseq)])+'\n')
    def test_softmask_preserves_bytes(self):
        p = self.file('g.fa', '>lowercase_header\nACGT\n>two\nacgt\n')
        before = p.read_bytes(); self.assertTrue(guards.softmask(p)); self.assertEqual(p.read_bytes(), before)
    def test_lowercase_header_is_not_masking(self):
        p = self.file('g.fa', '>lowercase\nACGTN\n')
        with self.assertRaisesRegex(ValueError, 'No lowercase'): guards.softmask(p)
    def test_chunk_boundary_header(self):
        p = self.file('g.fa', '>' + 'a' * (1024**2+5) + '\nACGT\n')
        with self.assertRaisesRegex(ValueError, 'No lowercase'): guards.softmask(p)
    def test_ep_manifest(self):
        s = guards.samples(self.manifest())[0]; self.assertEqual(s['runs'], [])
    def test_etp_multiple_pairs(self):
        rows = ['run_id\tfastq_r1\tfastq_r2']
        for i in range(2):
            rows.append('\t'.join([f'run{i}', str(self.file(f'{i}_1.fq')), str(self.file(f'{i}_2.fq'))]))
        rna = self.file('rna.tsv', '\n'.join(rows)+'\n')
        self.assertEqual(len(guards.samples(self.manifest('etp', rna))[0]['runs']), 2)
    def test_ep_rejects_rna(self):
        with self.assertRaises(ValueError): guards.samples(self.manifest('ep', self.file('rna.tsv')))
    def test_etp_requires_distinct_mates(self):
        fastq = self.file('reads.fq')
        rna = self.file('rna.tsv', f'run_id\tfastq_r1\tfastq_r2\nr1\t{fastq}\t{fastq}\n')
        with self.assertRaises(ValueError): guards.samples(self.manifest('etp', rna))
    def test_duplicate_sample(self):
        p = self.manifest(); text = p.read_text(); p.write_text(text + text.splitlines()[1]+'\n')
        with self.assertRaises(ValueError): guards.samples(p)
    def test_fingerprint_change(self):
        s = guards.samples(self.manifest())[0]; old = guards.fingerprints(s)
        Path(s['protein_fasta']).write_text('>p1\nMBBBB\n')
        self.assertNotEqual(old, guards.fingerprints(s))
    def test_reject_unowned_root(self):
        out, work = self.root/'out', self.root/'scratch'
        out.mkdir(); (out/'important').write_text('keep')
        with self.assertRaises(ValueError): guards.claim_roots(out, work)
        self.assertTrue((out/'important').exists())
    def test_root_binding_and_idempotence(self):
        out, work = self.root/'out', self.root/'scratch'
        guards.claim_roots(out, work); guards.claim_roots(out, work)
        with self.assertRaises(ValueError): guards.claim_roots(out, self.root/'different')
    def test_symlink_owned_root(self):
        target = self.root/'target'; target.mkdir(); link = self.root/'link'; link.symlink_to(target)
        with self.assertRaises(ValueError): guards.workspace(str(link))
    def test_space_check_readonly(self):
        cfg = dict(work_root=str(self.root/'scratch'), output_root=str(self.root/'out'),
                   min_free_work_gb=0, min_free_output_gb=0)
        self.assertEqual(set(driver.resources(cfg)), {'work_root', 'output_root'})
        self.assertFalse((self.root/'scratch').exists())
        cfg['min_free_work_gb'] = 10**9
        with self.assertRaises(ValueError): driver.resources(cfg)


if __name__ == '__main__': unittest.main()
