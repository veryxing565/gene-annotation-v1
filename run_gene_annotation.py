#!/usr/bin/env python3
"""Serial arachnid protein-coding annotation: fastp, BRAKER4, nine filters, compleasm."""
from __future__ import annotations
import argparse
import configparser
import contextlib
import csv
from datetime import datetime
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE / 'scripts'))
import runtime as guards


class Cancelled(InterruptedError):
    def __init__(self, signum):
        self.exit_code = 128 + signum
        super().__init__('Received signal ' + str(signum))


def now(): return datetime.now().astimezone().isoformat()


def core():
    sys.path.insert(0, str(BASE / 'pipeline'))
    import pipeline
    return pipeline


def configuration(path):
    cfg = json.loads(Path(path).read_text(encoding='utf-8'))
    cfg['output_root'] = str(guards.workspace(cfg['output_root']))
    cfg['work_root'] = str(guards.workspace(cfg['work_root']))
    output, work = Path(cfg['output_root']), Path(cfg['work_root'])
    guards.require(output != work and output not in work.parents and work not in output.parents, 'Output/work roots must be disjoint')
    guards.require(work != BASE.resolve() and output != BASE.resolve() and work not in BASE.resolve().parents and output not in BASE.resolve().parents and BASE.resolve() not in work.parents and BASE.resolve() not in output.parents, 'Keep work/output separate from the installed source')
    cfg.setdefault('threads', 56); cfg.setdefault('memory_mb', 240000)
    cfg.setdefault('rule_memory_mb', min(120000, cfg['memory_mb']))
    cfg.setdefault('fastp_threads', min(8, cfg['threads']))
    cfg.setdefault('cleanup_work', True); cfg.setdefault('max_consecutive_failures', 3)
    cfg.setdefault('min_free_work_gb', 100); cfg.setdefault('min_free_output_gb', 20)
    for key in ['threads', 'fastp_threads', 'memory_mb', 'rule_memory_mb', 'max_consecutive_failures']:
        guards.require(type(cfg[key]) is int and cfg[key] > 0, 'Positive integer required: ' + key)
    for key in ['min_free_work_gb', 'min_free_output_gb']:
        guards.require(type(cfg[key]) in (int, float) and 0 <= cfg[key] < float('inf'), 'Finite nonnegative space limit required: ' + key)
    guards.require(1 <= cfg['threads'] <= 56 and cfg['fastp_threads'] <= cfg['threads'], 'CPU budget must be 1..56; fastp cannot exceed it')
    guards.require(cfg['rule_memory_mb'] <= cfg['memory_mb'], 'Rule memory must not exceed scheduling budget')
    guards.require(type(cfg['cleanup_work']) is bool, 'cleanup_work must be boolean')
    cfg['tools'] = dict(cfg['tools'])
    for key in ['snakemake', 'apptainer', 'compleasm', 'fastp']:
        value = cfg['tools'].get(key, key)
        candidate = value if '/' in value else shutil.which(value)
        guards.require(bool(candidate), 'Missing executable: ' + key)
        p = guards.safe_path(candidate, existing=True)
        guards.require(os.access(p, os.X_OK), 'Executable permission missing: ' + str(p))
        cfg['tools'][key] = str(p)
    for key in ['braker3_image', 'agat_image', 'busco_image']:
        cfg['containers'][key] = str(guards.safe_path(cfg['containers'][key], existing=True))
    cfg['databases']['compleasm'] = str(guards.safe_path(cfg['databases']['compleasm'], existing=True, directory=True))
    cfg['databases']['busco'] = str(guards.safe_path(cfg['databases'].get('busco', cfg['databases']['compleasm']), existing=True, directory=True))
    for key in ['container_cache', 'apptainer_cache']:
        cfg[key] = str(guards.safe_path(cfg[key]))
        p = Path(cfg[key])
        guards.require(p != work and work not in p.parents, 'Dependency cache must not be inside scratch')
    cfg.setdefault('java_home', '')
    if cfg['java_home']: cfg['java_home'] = str(guards.safe_path(cfg['java_home'], existing=True, directory=True))
    for value in [*cfg['tools'].values(), *cfg['containers'].values(), *cfg['databases'].values()]:
        p = Path(value)
        guards.require(p != work and work not in p.parents, 'Keep dependencies outside the cleanup workspace')
    return cfg


def resources(cfg):
    """Free-space thresholds are preflight checks, not space reservations."""
    free = {}
    for location, limit in [('work_root', 'min_free_work_gb'), ('output_root', 'min_free_output_gb')]:
        parent = Path(cfg[location])
        while not parent.exists(): parent = parent.parent
        available = shutil.disk_usage(parent).free / 1024**3
        guards.require(available >= cfg[limit], 'Insufficient free space for ' + location + ': {:.1f} GiB'.format(available))
        free[location] = available
    return free


def check(cfg, samples, versions=True):
    guards.require(sys.platform.startswith('linux') and not sys.flags.optimize, 'Linux and Python without -O/PYTHONOPTIMIZE are required')
    guards.safe_path(str(BASE.resolve()), existing=True, directory=True)
    guards.require(cfg['threads'] <= len(os.sched_getaffinity(0)), 'Requested CPU count exceeds scheduler affinity')
    resources(cfg)
    work = Path(cfg['work_root'])
    output = Path(cfg['output_root'])
    for sample in samples:
        for p in guards.input_paths(sample):
            guards.require(work != p and work not in p.parents and output != p and output not in p.parents, 'Keep inputs outside dedicated work/output roots')
        guards.softmask(sample['genome'])
        core().assessment(Path(sample['genome_compleasm_summary']).parent, 'arachnida')
    for lineage in ['arachnida', 'arthropoda']:
        db = Path(cfg['databases']['compleasm']) / (lineage + '_odb12.2')
        for name in ['lengths_cutoff', 'scores_cutoff', 'refseq_db.faa.gz']:
            guards.require((db / name).is_file(), 'Missing compleasm cache: ' + str(db / name))
        guards.require((db / 'hmms').is_dir() and any((db / 'hmms').iterdir()), 'Missing HMM cache')
    found = {}
    if versions:
        for key, expected in [('snakemake', '8.18.2'), ('compleasm', '0.2.9'), ('apptainer', None), ('fastp', '0.23.3')]:
            r = subprocess.run([cfg['tools'][key], '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=30)
            guards.require(r.returncode == 0 and (expected is None or expected in r.stdout), 'Wrong/unavailable ' + key + ': ' + r.stdout)
            found[key] = r.stdout.strip()
    return found


def source_hashes():
    paths = [BASE / 'run_gene_annotation.py', BASE / 'scripts/runtime.py', BASE / 'config/braker.parameters.ini']
    paths += list((BASE / 'pipeline').glob('*.py'))
    paths += [p for p in (BASE / 'workflow').rglob('*') if p.is_file() and
              (p.name == 'Snakefile' or p.suffix in ('.smk', '.py', '.pl', '.sh'))]
    return {str(p.relative_to(BASE)): guards.sha(p) for p in sorted(paths)}


def native_call(argv, log, environment):
    fds = tuple(int(x) for x in environment['GENEANNOT_LEASE_FDS'].split(','))
    with Path(log).open('a', encoding='utf-8') as stream:
        stream.write(now() + ' COMMAND ' + repr(argv) + '\n'); stream.flush()
        child = subprocess.Popen([str(x) for x in argv], stdin=subprocess.DEVNULL, stdout=stream,
                                 stderr=subprocess.STDOUT, env=environment, start_new_session=True, pass_fds=fds)
        try: code = child.wait()
        except BaseException:
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGTERM)
                try: child.wait(timeout=15)
                except subprocess.TimeoutExpired: os.killpg(child.pid, signal.SIGKILL); child.wait()
            raise
        guards.require(code == 0, 'External tool failed; inspect ' + str(log))


def quality_control(cfg, sample, output, temporary, environment, update):
    cleaned = dict(sample=sample['sample_id'], r1=[], r2=[], records=[])
    qwork = temporary / 'qc'; qwork.mkdir()
    archived = output / '00_config/rna_qc'; archived.mkdir(exist_ok=True)
    for run in sample['runs']:
        rid = run['run_id']; update('fastp', run_id=rid)
        r1 = qwork / (rid + '.R1.fastq.gz'); r2 = qwork / (rid + '.R2.fastq.gz')
        js = qwork / (rid + '.json'); html = qwork / (rid + '.html')
        native_call([cfg['tools']['fastp'], '--in1', run['fastq_r1'], '--in2', run['fastq_r2'],
                     '--out1', r1, '--out2', r2, '--json', js, '--html', html,
                     '--thread', str(cfg['fastp_threads']), *guards.FASTP_ARGS], output / 'logs' / (rid + '.fastp.log'), environment)
        report = json.loads(js.read_text())
        guards.require(report['summary']['after_filtering']['total_reads'] > 0, 'No reads survived fastp')
        for p in [r1, r2]: guards.require(p.is_file() and p.stat().st_size > 0, 'Missing clean FASTQ')
        cleaned['r1'].append(str(r1)); cleaned['r2'].append(str(r2))
        for p in [js, html]: core().copy(p, archived / p.name)
        cleaned['records'].append(dict(run_id=rid, summary=report['summary'],
                                       clean_sha256={str(r1):guards.sha(r1), str(r2):guards.sha(r2)}))
    guards.write_json(archived / 'summary.json', cleaned)
    return cleaned


def sample_config(cfg, sample, output, temporary, cleaned):
    configdir = output / '00_config'
    template = configparser.ConfigParser(); template.read(BASE / 'config/braker.parameters.ini')
    template['paths'] = dict(samples_file=str(configdir / 'samples.csv'),
                             augustus_config_path=str(temporary / '01_braker/augustus_config'),
                             busco_download_path=cfg['databases']['busco'], compleasm_download_path=cfg['databases']['compleasm'])
    template['containers'] = cfg['containers']
    template['SLURM_ARGS'] = dict(cpus_per_task=str(cfg['threads']), mem_of_node=str(cfg['rule_memory_mb']), max_runtime='4320')
    ini = configdir / 'braker.template.ini'
    with ini.open('w') as f: template.write(f)
    fields = 'sample_name genome genome_masked protein_fasta bam_files fastq_r1 fastq_r2 sra_ids varus_genus varus_species isoseq_bam isoseq_fastq busco_lineage'.split()
    row = dict.fromkeys(fields, '')
    row.update(sample_name=sample['sample_id'], genome=sample['genome'], genome_masked=sample['genome'],
               protein_fasta=sample['protein_fasta'], fastq_r1=':'.join(cleaned['r1']), fastq_r2=':'.join(cleaned['r2']), busco_lineage='arthropoda_odb12.2')
    csv_path = configdir / 'samples.csv'
    with csv_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader(); writer.writerow(row)
    bind = {BASE.resolve(), output, temporary, *[p.parent for p in guards.input_paths(sample)],
            Path(cfg['databases']['compleasm']), Path(cfg['databases']['busco'])}
    runtime = {k:cfg['tools'][k] for k in ['snakemake', 'apptainer', 'compleasm']}
    runtime.update(compleasm_db=cfg['databases']['compleasm'], container_cache=cfg['container_cache'],
                   apptainer_cache=cfg['apptainer_cache'], apptainer_tmpdir=str(temporary / 'container_tmp'),
                   java_home=cfg['java_home'], memory_mb=cfg['memory_mb'], bind_paths=sorted(map(str, bind)))
    return dict(sample=sample['sample_id'], mode=sample['mode'], output_root=str(output), temporary_root=str(temporary),
                threads=cfg['threads'], braker_template=str(ini), samples_csv=str(csv_path),
                workflow=str(BASE / 'workflow/Snakefile'), genome_compleasm_summary=sample['genome_compleasm_summary'], runtime=runtime)


def run_sample(cfg, sample, batch_fd, versions):
    import fcntl
    output = guards.sample_output(Path(cfg['output_root']), sample['sample_id'])
    temporary = Path(cfg['work_root']) / sample['sample_id']
    lock_path = output / 'logs/driver.lock'; guards.no_links(lock_path)
    with lock_path.open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        state = dict(sample=sample['sample_id'], mode=sample['mode'], started=now(), threads=cfg['threads'],
                     cpu_affinity=sorted(os.sched_getaffinity(0)))
        def update(stage, **extra):
            state.update(stage=stage, updated=now(), **extra)
            guards.write_json(output / 'logs/driver.json', state)
            print(now(), sample['sample_id'], stage, flush=True)
        update('fingerprint_inputs')
        inputs = guards.fingerprints(sample)
        sources = source_hashes()
        provenance = dict(sample=sample, configuration=cfg, inputs=inputs, versions=versions, workflow_source_sha256=sources)
        existing = output / '00_config/input_provenance.json'
        if (output / 'pipeline.complete.json').exists():
            old = json.loads(existing.read_text())
            guards.require(old['inputs'] == inputs and old['configuration'] == cfg and old['sample'] == sample and old['workflow_source_sha256'] == sources, 'Completed sample inputs/config/source changed')
            core().verify_completed(output)
            update('skipped_verified_complete'); return 'skipped'
        guards.require(not existing.exists() and not temporary.exists(), 'Incomplete sample exists; inspect it before a fresh retry')
        resources(cfg)
        guards.no_links(temporary); temporary.mkdir()
        (temporary / '.gene-driver-owned').write_text(str(output) + '\n')
        (temporary / 'container_tmp').mkdir()
        guards.write_json(existing, provenance)
        for key in ['container_cache', 'apptainer_cache']: Path(cfg[key]).mkdir(parents=True, exist_ok=True)
        environment = {k:v for k,v in os.environ.items() if not k.startswith('BRAKER4_') and k != 'GENEANNOT_LEASE_FDS'}
        environment.update(GENEANNOT_LEASE_FDS=','.join(map(str, [batch_fd, lock.fileno()])),
                           OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONNOUSERSITE='1')
        previous_lease = os.environ.get('GENEANNOT_LEASE_FDS')
        try:
            os.environ['GENEANNOT_LEASE_FDS'] = environment['GENEANNOT_LEASE_FDS']
            cleaned = quality_control(cfg, sample, output, temporary, environment, update) if sample['mode'] == 'etp' else dict(r1=[], r2=[])
            pipeline_cfg = sample_config(cfg, sample, output, temporary, cleaned)
            update('BRAKER4_nine_filters_compleasm')
            with (output / 'logs/pipeline.log').open('a') as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                core().run(pipeline_cfg)
            core().verify_completed(output)
            guards.require(guards.fingerprints(sample) == inputs, 'Input changed while annotation was running')
            guards.require(source_hashes() == sources, 'Workflow source changed while annotation was running')
            if cfg['cleanup_work']:
                guards.cleanup(Path(cfg['work_root']), output, temporary, core().verify_completed)
                guards.write_json(output / 'temporary_cleanup.json', dict(removed=str(temporary), time=now(), checksums_verified=True))
            update('completed', temporary_removed=cfg['cleanup_work'])
            return 'completed'
        except BaseException as exc:
            update('failed_temporary_preserved', error=repr(exc))
            raise
        finally:
            if previous_lease is None: os.environ.pop('GENEANNOT_LEASE_FDS', None)
            else: os.environ['GENEANNOT_LEASE_FDS'] = previous_lease


def run(cfg, samples, versions):
    import fcntl
    output, work = Path(cfg['output_root']), Path(cfg['work_root'])
    guards.claim_roots(output, work)
    lock_path = output / 'state/batch.lock'; guards.no_links(lock_path)
    with lock_path.open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        cpus = sorted(os.sched_getaffinity(0))[:cfg['threads']]; os.sched_setaffinity(0, cpus)
        state = dict(started=now(), pid=os.getpid(), cpu_affinity=cpus, threads=cfg['threads'], complete=[], skipped=[], failed=[])
        consecutive = 0
        def interrupted(signum, frame):
            state.update(stage='interrupted_files_preserved', finished=now(), returncode=128 + signum)
            guards.write_json(output / 'state/batch.json', state)
            raise Cancelled(signum)
        signal.signal(signal.SIGTERM, interrupted); signal.signal(signal.SIGINT, interrupted)
        for sample in samples:
            state.update(current=sample['sample_id'], stage='running', updated=now())
            guards.write_json(output / 'state/batch.json', state)
            try:
                result = run_sample(cfg, sample, lock.fileno(), versions)
                state[result if result == 'skipped' else 'complete'].append(sample['sample_id'])
                if result != 'skipped': consecutive = 0
            except (KeyboardInterrupt, InterruptedError): raise
            except Exception as exc:
                state['failed'].append(dict(sample=sample['sample_id'], error=repr(exc)))
                consecutive += 1
                print('FAILED', sample['sample_id'], repr(exc), file=sys.stderr, flush=True)
            if consecutive >= cfg['max_consecutive_failures']: break
        code = 3 if consecutive >= cfg['max_consecutive_failures'] else 1 if state['failed'] else 0
        state.update(stage='completed' if code == 0 else 'finished_with_failures', current=None, finished=now(), returncode=code)
        guards.write_json(output / 'state/batch.json', state)
        if code == 0 and cfg['cleanup_work'] and set(work.iterdir()) == {work / '.gene-v1-root.json'}:
            (work / '.gene-v1-root.json').unlink(); work.rmdir()
        return code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=(BASE / 'VERSION').read_text().strip())
    parser.add_argument('action', choices=['check', 'run'])
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--samples', type=Path, required=True)
    args = parser.parse_args()
    cfg, selected = configuration(args.config), guards.samples(args.samples)
    versions = check(cfg, selected)
    if args.action == 'check':
        print(json.dumps(dict(ready=True, samples=[s['sample_id'] for s in selected], versions=versions), indent=2)); return 0
    return run(cfg, selected, versions)


if __name__ == '__main__':
    try: raise SystemExit(main())
    except Cancelled as exc:
        print('INTERRUPTED:', exc, file=sys.stderr); raise SystemExit(exc.exit_code)
    except (ValueError, OSError, AssertionError) as exc:
        print('ERROR:', exc, file=sys.stderr); raise SystemExit(2)
