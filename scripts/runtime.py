"""Validation, ownership and input fingerprints; no biological prediction."""
from __future__ import annotations
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil

SCHEMA = 'gene-annotation-v1-portable-1'
FASTP_ARGS = ['--detect_adapter_for_pe', '--trim_poly_g', '--qualified_quality_phred', '20',
              '--unqualified_percent_limit', '40', '--n_base_limit', '5', '--length_required', '50']
SAMPLE_COLUMNS = ['sample_id', 'mode', 'genome', 'protein_fasta', 'genome_compleasm_summary', 'rnaseq_manifest']


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''): h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path); partial = path.with_name(path.name + '.partial')
    no_links(path); no_links(partial)
    partial.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    partial.replace(path)


def no_links(path):
    path = Path(path)
    require(path.is_absolute(), 'Absolute path required: ' + str(path))
    for p in [path, *path.parents]: require(not p.is_symlink(), 'Owned path contains a symlink: ' + str(p))


def safe_path(value, existing=False, directory=False):
    path = Path(value)
    require(path.is_absolute() and '..' not in path.parts, 'Use an absolute path without ..: ' + str(path))
    require(bool(re.fullmatch(r'/[A-Za-z0-9_./+\-]+', str(path))), 'Shell-safe Linux paths required: ' + str(path))
    if existing:
        require(path.is_dir() if directory else path.is_file(), 'Missing path: ' + str(path))
        require(bool(re.fullmatch(r'/[A-Za-z0-9_./+\-]+', str(path.resolve()))), 'Resolved path is not shell-safe')
        return path.resolve()
    no_links(path)
    return path


def workspace(value):
    path = safe_path(value)
    broad = {Path(p) for p in ['/', '/home', '/tmp', '/var/tmp', '/mnt', '/media', '/scratch', '/data', '/work', '/opt']}
    broad.add(Path.home().resolve())
    require(path not in broad and len(path.parts) >= 3 and not path.is_mount(), 'Use a dedicated workspace subdirectory')
    return path


def table(path, columns):
    with Path(path).open(encoding='utf-8', newline='') as stream:
        lines = [line for line in stream if line.strip() and not line.startswith('#')]
    require(all(not line.endswith('\r\n') for line in lines), 'Use Unix LF TSV line endings')
    reader = csv.DictReader(lines, delimiter='\t')
    require(reader.fieldnames == columns, 'TSV columns must be exactly: ' + '\t'.join(columns))
    rows = list(reader)
    require(bool(rows), 'Empty manifest')
    require(all(set(row) == set(columns) and all(v is not None for v in row.values()) for row in rows), 'Malformed TSV row')
    return rows


def samples(path):
    result = []; seen = set()
    for row in table(path, SAMPLE_COLUMNS):
        name = row['sample_id']
        require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,119}', name)) and name not in seen, 'Unsafe/duplicate sample ID')
        seen.add(name)
        require(row['mode'] in ('ep', 'etp'), 'Only ep and etp are supported')
        for key in ['genome', 'protein_fasta', 'genome_compleasm_summary']:
            row[key] = str(safe_path(row[key], existing=True))
            require(Path(row[key]).stat().st_size > 0, 'Empty input: ' + row[key])
        require(Path(row['genome_compleasm_summary']).name == 'summary.txt', 'The genome assessment file must be named summary.txt')
        runs = []
        if row['mode'] == 'etp':
            manifest = safe_path(row['rnaseq_manifest'], existing=True)
            row['rnaseq_manifest'] = str(manifest)
            ids = set()
            for run in table(manifest, ['run_id', 'fastq_r1', 'fastq_r2']):
                require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,119}', run['run_id'])) and run['run_id'] not in ids, 'Unsafe/duplicate run ID')
                ids.add(run['run_id'])
                for key in ['fastq_r1', 'fastq_r2']:
                    run[key] = str(safe_path(run[key], existing=True))
                    require(Path(run[key]).stat().st_size > 0, 'Empty FASTQ')
                require(run['fastq_r1'] != run['fastq_r2'], 'R1 and R2 must be distinct files')
                runs.append(run)
        else:
            require(row['rnaseq_manifest'] in ('', '-'), 'EP must not silently ignore supplied RNA evidence')
        row['runs'] = runs
        result.append(row)
    return result


def input_paths(sample):
    paths = [Path(sample[k]) for k in ['genome', 'protein_fasta', 'genome_compleasm_summary']]
    if sample['runs']: paths.append(Path(sample['rnaseq_manifest']))
    paths += [Path(run[k]) for run in sample['runs'] for k in ['fastq_r1', 'fastq_r2']]
    return paths


def fingerprints(sample):
    return {str(p): dict(bytes=p.stat().st_size, sha256=sha(p)) for p in input_paths(sample)}


def softmask(path):
    """Bounded-memory case inspection, excluding headers even across chunks."""
    at_start = True; header = False; found_sequence = False
    with Path(path).open('rb') as stream:
        first = stream.read(1)
        require(first == b'>', 'Provide an uncompressed FASTA beginning with >')
        stream.seek(0)
        for chunk in iter(lambda: stream.read(1024**2), b''):
            for part in chunk.splitlines(keepends=True):
                if at_start: header = part.startswith(b'>')
                if not header:
                    found_sequence |= bool(re.search(rb'[ACGTNacgtn]', part))
                    if re.search(rb'[acgt]', part): return True
                at_start = part.endswith((b'\n', b'\r'))
                if at_start: header = False
    require(found_sequence, 'No DNA sequence found in FASTA')
    raise ValueError('No lowercase acgt found; provide a real soft-masked assembly')


def claim_roots(output, work):
    for path, kind in [(output, 'output'), (work, 'work')]:
        no_links(path)
        marker = path / '.gene-v1-root.json'
        expected = dict(schema=SCHEMA, kind=kind, output_root=str(output), work_root=str(work))
        if marker.exists():
            no_links(marker); require(json.loads(marker.read_text()) == expected, 'Workspace ownership mismatch')
        else:
            require(not path.exists() or not any(path.iterdir()), 'Refusing nonempty unowned workspace: ' + str(path))
            path.mkdir(parents=True, exist_ok=True); write_json(marker, expected)
    for name in ['results', 'state', 'logs']:
        p = output / name; no_links(p); p.mkdir(exist_ok=True)


def sample_output(output_root, name):
    output = output_root / 'results' / name; no_links(output)
    marker = output / '.gene-sample.json'
    expected = dict(schema=SCHEMA, sample=name, output=str(output))
    if marker.exists(): require(json.loads(marker.read_text()) == expected, 'Sample output ownership mismatch')
    else:
        require(not output.exists() or not any(output.iterdir()), 'Unowned existing sample output')
        output.mkdir(parents=True, exist_ok=True); write_json(marker, expected)
    for name in ['logs', '00_config']:
        no_links(output / name); (output / name).mkdir(exist_ok=True)
    return output


def cleanup(work, output, temporary, verify):
    no_links(temporary)
    require(temporary.parent == work and temporary.name == output.name, 'Unsafe cleanup target')
    require((temporary / '.gene-driver-owned').read_text().strip() == str(output), 'Driver ownership mismatch')
    require((temporary / '.pipeline-owned').read_text().strip() == str(output), 'Pipeline ownership mismatch')
    verify(output)
    require(shutil.rmtree.avoids_symlink_attacks, 'Safe fd-based cleanup is required')
    shutil.rmtree(temporary)
