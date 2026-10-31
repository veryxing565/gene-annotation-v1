#!/usr/bin/env python3
"""TSEBRA-only threshold trials using the unchanged BRAKER postprocessing rules."""
from __future__ import annotations

import configparser
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess

ROOT = TMP = BASE = CONFIG_DIR = WORKFLOW = None
SAMPLE = None
THREADS = 56
MODE = 'etp'
RUNTIME = {}
ENV = {}
APPT = []
ALLOWED = ['best_by_compleasm', 'both_agree_rescue', 'filter_internal_stop_codons',
           'normalize_cds', 'extract_final_sequences', 'add_utr', 'fix_gtf',
           'convert_gtf_to_gff3', 'gene_support_summary']
INPUTS = []


def now():
    return datetime.datetime.now().astimezone().isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def no_symlinks(path):
    path = Path(path)
    assert path.is_absolute(), 'Absolute path required: ' + str(path)
    assert path.resolve() == path, 'Noncanonical or symlink path: ' + str(path)
    for part in [path, *path.parents]:
        assert not part.is_symlink(), 'Symlink not permitted: ' + str(part)
    return path


def write_json(path, value):
    path = no_symlinks(Path(path))
    temporary = path.with_name(path.name + '.tmp')
    no_symlinks(temporary)
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    temporary.replace(path)


def read_json(path):
    return json.loads(no_symlinks(Path(path)).read_text(encoding='utf-8'))


def regular_tree(directory):
    directory = no_symlinks(Path(directory))
    assert directory.is_dir(), 'Directory missing: ' + str(directory)
    files = set()
    for path in directory.rglob('*'):
        assert not path.is_symlink(), 'Symlink in result tree: ' + str(path)
        assert path.is_dir() or path.is_file(), 'Nonregular result: ' + str(path)
        if path.is_file():
            no_symlinks(path)
            files.add(path)
    return files


def write_manifest(directory):
    directory = Path(directory)
    manifest = directory / 'SHA256SUMS'
    no_symlinks(manifest)
    files = regular_tree(directory) - {manifest}
    assert files, 'Cannot publish an empty checksum manifest'
    manifest.write_text(''.join(sha(path) + '  ' + path.relative_to(directory).as_posix() + '\n'
                                for path in sorted(files)), encoding='utf-8')
    verify_manifest(directory)


def verify_manifest(directory):
    """Verify exact SHA256 coverage; reject links, duplicates and path escapes."""
    directory = no_symlinks(Path(directory))
    manifest = no_symlinks(directory / 'SHA256SUMS')
    assert manifest.is_file() and manifest.stat().st_size, 'Missing/empty result SHA256SUMS'
    covered = set()
    for line in manifest.read_text(encoding='utf-8').splitlines():
        match = re.fullmatch(r'([0-9a-fA-F]{64}) [ *](.+)', line)
        assert match, 'Malformed checksum manifest line'
        expected, relative = match.groups()
        relative = Path(relative)
        assert not relative.is_absolute() and '..' not in relative.parts, 'Checksum path escape'
        path = no_symlinks(directory / relative)
        assert path != manifest and path not in covered, 'Self-reference or duplicate checksum entry'
        assert path.is_relative_to(directory) and path.is_file(), 'Missing checksum target'
        assert sha(path) == expected.lower(), 'Result checksum mismatch: ' + str(path)
        covered.add(path)
    assert covered, 'Empty result checksum manifest'
    assert covered == regular_tree(directory) - {manifest}, 'Checksum manifest does not cover all results'
    return len(covered)


def copy_file(source, target):
    """Materialize approved input file links; never write through destination links."""
    source, target = Path(source), Path(target)
    assert source.is_file(), 'Missing source file: ' + str(source)
    no_symlinks(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    expected = sha(source)
    if target.exists():
        assert target.is_file() and sha(target) == expected, 'Conflicting file: ' + str(target)
    else:
        shutil.copy2(source, target)
    assert sha(target) == expected, 'Copy checksum mismatch: ' + str(target)


def copy_tree(source, target):
    source, target = Path(source), no_symlinks(Path(target))
    files = regular_tree(source)
    if target.exists():
        assert target.is_dir()
        assert {p.relative_to(source) for p in files} == {p.relative_to(target) for p in regular_tree(target)}, 'Conflicting directory archive'
    else:
        target.mkdir(parents=True)
    for path in sorted(files):
        copy_file(path, target / path.relative_to(source))


def configure_runtime(runtime, braker_image):
    """Reset all process-local runtime configuration for each serialized sample."""
    global RUNTIME, ENV, APPT
    required = {'snakemake', 'apptainer', 'compleasm', 'compleasm_db', 'container_cache',
                'apptainer_cache', 'apptainer_tmpdir', 'memory_mb', 'bind_paths'}
    assert required <= set(runtime), 'Missing runtime keys: ' + repr(required - set(runtime))
    current = dict(runtime)
    for key in ('snakemake', 'apptainer', 'compleasm'):
        path = Path(current[key])
        assert path.is_absolute() and path.is_file() and os.access(path, os.X_OK), 'Invalid executable: ' + key
        current[key] = str(path)
    database = Path(current['compleasm_db'])
    assert database.is_absolute() and database.is_dir(), 'Missing local compleasm library directory'
    current['compleasm_db'] = str(database)
    for key in ('container_cache', 'apptainer_cache', 'apptainer_tmpdir'):
        path = no_symlinks(Path(current[key]))
        assert path not in {Path('/'), Path.home().resolve()}, 'Use a dedicated cache directory'
        path.mkdir(parents=True, exist_ok=True)
        current[key] = str(path)
    current['memory_mb'] = int(current['memory_mb'])
    assert current['memory_mb'] > 0, 'memory_mb must be positive'
    bindings = []
    assert isinstance(current['bind_paths'], (list, tuple)) and current['bind_paths'], 'Explicit bind_paths required'
    for raw in current['bind_paths']:
        path = Path(raw)
        assert path.is_absolute() and path.is_dir(), 'Bind directory must already exist: ' + str(path)
        assert not any(c in str(path) for c in ',:\r\n'), 'Unsupported Apptainer bind delimiter'
        if str(path) not in bindings:
            bindings.append(str(path))
    current['bind_paths'] = bindings
    current['bind_argument'] = ','.join(bindings)
    current['apptainer_args'] = shlex.join(['-B', current['bind_argument']])
    image = Path(braker_image)
    assert image.is_absolute() and image.is_file(), 'A local BRAKER image is required'
    current['braker_image'] = str(image)
    env = {k: v for k, v in os.environ.items() if not k.startswith('BRAKER4_')}
    tool_directories = list(dict.fromkeys(str(Path(current[key]).parent) for key in ('snakemake', 'apptainer', 'compleasm')))
    env.update(PATH=':'.join(tool_directories + [env.get('PATH', '')]), PYTHONNOUSERSITE='1',
               APPTAINER_CACHEDIR=current['apptainer_cache'], APPTAINER_TMPDIR=current['apptainer_tmpdir'],
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    java_home = current.get('java_home')
    if java_home:
        assert Path(java_home).is_absolute() and Path(java_home).is_dir(), 'Invalid java_home'
        env['JAVA_HOME'] = str(java_home)
    # All supported datasets/images are local. Never silently use a user proxy
    # to retrieve missing data while a scientific run is in progress.
    for key in ('http_proxy', 'https_proxy', 'HTTP_PROXY', 'HTTPS_PROXY', 'all_proxy', 'ALL_PROXY'):
        env[key] = 'http://127.0.0.1:9'
    env.update(no_proxy='', NO_PROXY='')
    RUNTIME, ENV = current, env
    APPT = [current['apptainer'], 'exec', '--bind', current['bind_argument'], str(image)]


def status(trial, stage):
    write_json(ROOT.parent / 'logs/current.json', dict(trial=trial, stage=stage, updated=now()))
    print(now(), trial, stage, flush=True)


def lease_fds():
    raw = os.environ.get('GENEANNOT_LEASE_FDS', '')
    if not raw:
        return ()
    values = raw.split(',')
    assert all(v.isdecimal() for v in values), 'Invalid GENEANNOT_LEASE_FDS'
    result = tuple(sorted(set(map(int, values))))
    for fd in result:
        assert fd >= 3, 'Only explicitly opened lease descriptors may be inherited'
        os.fstat(fd)
    return result


def call(command, log, env=None, cwd=None):
    assert RUNTIME, 'configure_runtime() must precede scientific commands'
    log = no_symlinks(Path(log))
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open('w', encoding='utf-8') as stream:
        stream.write('COMMAND ' + repr(command) + '\n')
        stream.flush()
        child = subprocess.Popen([str(x) for x in command], env=ENV if env is None else env,
                                 cwd=cwd, stdin=subprocess.DEVNULL, stdout=stream,
                                 stderr=subprocess.STDOUT, start_new_session=True,
                                 pass_fds=lease_fds())
        try:
            code = child.wait()
        except BaseException:
            if child.poll() is None:
                try:
                    os.killpg(child.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    child.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(child.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    child.wait(timeout=15)
            raise
        if code:
            raise subprocess.CalledProcessError(code, command)


def snakemake_runtime_arguments():
    return ['--resources', 'mem_mb=' + str(RUNTIME['memory_mb']), '--use-apptainer',
            '--apptainer-prefix', RUNTIME['container_cache'],
            '--apptainer-args', RUNTIME['apptainer_args'], '--latency-wait', '120', '--restart-times', '0']


def trial_owner(work, label):
    return dict(schema='gene-annotation-trial-v1', sample=SAMPLE, trial=label,
                output_root=str(ROOT.parent), temporary_root=str(TMP.parent), work=str(work))


def verify_case(case, label, intron, stasto, single):
    case = no_symlinks(Path(case))
    assert (case / 'complete.json').is_file(), 'Trial completion marker missing'
    completion = read_json(case / 'complete.json')
    assert completion['results'] == str(case / 'results'), 'Trial result directory changed'
    parameters = read_json(case / 'parameters.json')
    assert (parameters['intron_support'], parameters['stasto_support'], parameters['filter_single_exon_genes']) == (float(intron), float(stasto), single), 'Trial parameters changed'
    checks = read_json(case / 'upstream_sha256.json')
    assert checks == {rel: sha(BASE / rel) for rel in INPUTS}, 'Different upstream evidence'
    verify_manifest(case / 'results')


def one(label, intron, stasto, single):
    assert label in 'ABCDEFGH' and len(label) == 1
    case = ROOT / label
    if (case / 'complete.json').exists():
        verify_case(case, label, intron, stasto, single)
        return
    no_symlinks(case)
    case.mkdir(exist_ok=True)
    logs = case / 'logs'
    no_symlinks(logs)
    logs.mkdir(exist_ok=True)
    work = no_symlinks(TMP / label)
    assert not work.exists(), 'Existing trial workspace: inspect before retry'
    assert (TMP.parent / '.pipeline-owned').read_text().strip() == str(ROOT.parent), 'Temporary pipeline ownership mismatch'
    dest = work / 'output' / SAMPLE
    dest.mkdir(parents=True)
    write_json(work / '.filter-trial-owned', trial_owner(work, label))
    status(label, 'stage_identical_upstream_inputs')
    checks = {}
    for relative in INPUTS:
        source, target = BASE / relative, dest / relative
        copy_file(source, target)
        checks[relative] = sha(source)
    write_json(case / 'upstream_sha256.json', checks)
    cfg = configparser.ConfigParser()
    cfg.read(str(CONFIG_DIR / 'config.ini'))
    assert cfg.getint('PARAMS', 'run_best_by_compleasm') == 0
    cfg['paths']['augustus_config_path'] = str(work / 'augustus_config')
    with (case / 'braker.config.ini').open('w') as stream:
        cfg.write(stream)
    tsebra_cfg = case / 'tsebra.cfg'
    tsebra_cfg.write_text('P 1\nE 20\nC 1\nM 1\nintron_support ' + intron +
                          '\nstasto_support ' + stasto + '\ne_1 0.1\ne_2 0.5\ne_3 0.05\ne_4 0.2\n')
    write_json(case / 'parameters.json', dict(label=label, intron_support=float(intron),
               stasto_support=float(stasto), filter_single_exon_genes=single, threads=THREADS,
               run_best_by_compleasm=False, workflow_source=str(WORKFLOW), work=str(work),
               postprocessing_rules=ALLOWED))
    status(label, 'TSEBRA')
    unrenamed = dest / 'tsebra.unrenamed.gtf'
    gm = 'GeneMark-EP/genemark.gtf' if MODE == 'ep' else 'GeneMark-ETP/genemark.gtf'
    keep = 'genemark/genemark.f.good.gtf' if MODE == 'ep' else 'GeneMark-ETP/training.gtf'
    hints = 'hintsfile_iter2.gff' if MODE == 'ep' else 'hintsfile.gff'
    if MODE == 'ep':
        # Preserve the production correction: raw and filtered GeneMark IDs can
        # otherwise collide after independent sorting/renumbering.
        for source, prefix, target in [(gm, 'EPRAW_', 'ep.raw.namespaced.gtf'),
                                       (keep, 'EPHC_', 'ep.good.namespaced.gtf')]:
            with (dest / source).open() as inp, (dest / target).open('w') as out:
                for line in inp:
                    out.write(re.sub(r'((?:gene_id|transcript_id) ")([^"]+)(")',
                                     lambda m: m[1] + prefix + m[2] + m[3], line))
        gm, keep = 'ep.raw.namespaced.gtf', 'ep.good.namespaced.gtf'
    command = APPT + ['tsebra.py', '--gtf', str(dest / 'augustus.hints.fixed.gtf') + ',' + str(dest / gm),
                      '--keep_gtf', str(dest / keep), '--hintfiles', str(dest / hints),
                      '--cfg', str(tsebra_cfg), '--out', str(unrenamed)]
    if single:
        command += ['--filter_single_exon_genes']
    call(command, logs / 'tsebra.log', cwd=str(work))
    call(APPT + ['/opt/TSEBRA/bin/rename_gtf.py', '--gtf', str(unrenamed),
                 '--out', str(dest / 'braker.tsebra.raw.gtf')], logs / 'rename.log', cwd=str(work))
    assert (dest / 'braker.tsebra.raw.gtf').stat().st_size > 0
    status(label, 'postprocessing_preflight')
    env = ENV.copy()
    env['BRAKER4_CONFIG'] = str(case / 'braker.config.ini')
    targets = [('output/' + SAMPLE + '/') + name for name in ('braker.aa', 'braker.gff3', 'gene_support.tsv')]
    command = [RUNTIME['snakemake']] + targets + ['--snakefile', str(WORKFLOW), '--directory', str(work),
               '--cores', str(THREADS)] + snakemake_runtime_arguments() + [
               '--rerun-triggers', 'mtime', '--printshellcmds', '--allowed-rules'] + ALLOWED
    call(command + ['--dry-run'], logs / 'postprocessing.dryrun.log', env)
    status(label, 'original_BRAKER_postprocessing')
    call(command, logs / 'postprocessing.log', env)
    assert (dest / 'braker.aa').stat().st_size > 0
    for lineage in ('arachnida', 'arthropoda'):
        status(label, 'compleasm_' + lineage)
        out = work / ('compleasm_' + lineage)
        call([RUNTIME['compleasm'], 'protein', '-p', str(dest / 'braker.aa'), '-l', lineage,
              '--odb', 'odb12.2', '-t', str(THREADS), '-L', RUNTIME['compleasm_db'], '-o', str(out)],
             logs / ('compleasm_' + lineage + '.log'))
        assert (out / 'summary.txt').is_file()
    status(label, 'transfer_and_validate')
    partial, final = case / 'results.partial', case / 'results'
    assert not partial.exists() and not partial.is_symlink() and not final.exists() and not final.is_symlink()
    partial.mkdir()
    (partial / 'annotation').mkdir()
    for name in ('braker.tsebra.raw.gtf', 'braker.tsebra.gtf', 'braker.rescued.gtf',
                 'braker.filtered.gtf', 'braker.gtf', 'braker.aa', 'braker.codingseq',
                 'braker_utr.gtf', 'braker.gff3', 'gene_support.tsv', 'normalize_cds.log',
                 'filter_stop_codons.log', 'both_agree_rescue.log', 'best_by_compleasm.log'):
        if (dest / name).is_file():
            copy_file(dest / name, partial / 'annotation' / name)
    for lineage in ('arachnida', 'arthropoda'):
        copy_tree(work / ('compleasm_' + lineage), partial / ('compleasm_' + lineage))
    write_manifest(partial)
    partial.rename(final)
    verify_manifest(final)
    for relative, expected in checks.items():
        assert sha(BASE / relative) == expected, 'Shared baseline was altered'
    write_json(case / 'complete.json', dict(completed=now(), results=str(final), work_retained=str(work)))
    status(label, 'complete')


def set_mode(mode):
    global MODE, INPUTS
    assert mode in ('ep', 'etp')
    MODE = mode
    INPUTS = ['genome.fa', 'genome.fa.fai', 'augustus.hints.fixed.gtf', 'hintsfile.gff']
    INPUTS += (['hintsfile_iter2.gff', 'GeneMark-EP/genemark.gtf', 'genemark/genemark.f.good.gtf']
               if mode == 'ep' else ['GeneMark-ETP/genemark.gtf', 'GeneMark-ETP/training.gtf',
                                    'GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff'])


def configure(root, temporary, baseline, config_dir, workflow, sample, threads):
    global ROOT, TMP, BASE, CONFIG_DIR, WORKFLOW, SAMPLE, THREADS
    ROOT, TMP, BASE, CONFIG_DIR = map(Path, (root, temporary, baseline, config_dir))
    WORKFLOW, SAMPLE, THREADS = Path(workflow), sample, int(threads)
    assert 1 <= THREADS <= 56
    no_symlinks(ROOT)
    no_symlinks(TMP)
    assert ROOT.parent != TMP.parent and TMP.parent not in ROOT.parents and ROOT.parent not in TMP.parents
    ROOT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
