#!/usr/bin/env python3
"""Sequential BRAKER4, nine TSEBRA settings, text reports and verified selection."""
from __future__ import annotations

import argparse
from collections import Counter
import configparser
import csv
import datetime
import json
import os
from pathlib import Path
import re
import shutil
import statistics
import sys

try:
    from . import filter_engine as engine
except ImportError:
    import filter_engine as engine

VERSION = 'gene-annotation-v1.0.0rc1'
PANEL = [('baseline', 1.0, 2, True), ('A', .9, 2, True), ('B', .8, 2, True),
         ('C', .8, 1, True), ('D', .8, 1, False), ('E', .75, 2, True),
         ('F', .70, 2, True), ('G', .60, 2, True), ('H', .50, 2, True)]
KEEP = ['braker.gtf', 'braker.aa', 'braker.codingseq', 'braker.gff3', 'braker_utr.gtf',
        'gene_support.tsv', 'both_agree_rescue.log', 'filter_stop_codons.log',
        'normalize_cds.log', 'best_by_compleasm.log']
REQUIRED_ANNOTATION = ['braker.gtf', 'braker.aa', 'braker.codingseq', 'braker.gff3', 'gene_support.tsv']
write_json = engine.write_json
read_json = engine.read_json
copy = engine.copy_file
verify_manifest = engine.verify_manifest


def stamp():
    return datetime.datetime.now().astimezone().isoformat()


def lines(path):
    with Path(path).open(encoding='utf-8') as stream:
        yield from stream


def assessment(path, lineage):
    path = Path(path)
    text = (path / 'summary.txt').read_text()
    found = re.search(r'lineage:\s*(\S+)', text)
    assert found and found.group(1) == lineage + '_odb12.2', str(path) + ' lineage mismatch'
    values = {key: int(value) for key, value in re.findall(r'^([SDFMI]):[^\n]*?,\s*(\d+)\s*$', text, re.M)}
    assert set(values) == set('SDFMI'), 'Missing compleasm categories'
    total = int(re.search(r'^N:\s*(\d+)', text, re.M).group(1))
    assert total > 0 and sum(values.values()) == total
    values.update(N=total, C=values['S'] + values['D'])
    values['C_percent'] = values['C'] * 100 / total
    table = path / 'full_table.tsv'
    if table.exists():
        identifiers = {}
        for line in lines(table):
            if line.startswith('#'):
                continue
            columns = line.rstrip().split('\t')
            if len(columns) >= 2:
                if columns[0] in identifiers:
                    assert identifiers[columns[0]] == columns[1]
                identifiers[columns[0]] = columns[1]
        assert len(identifiers) == total
        assert sum(value in ('Single', 'Duplicated') for value in identifiers.values()) == values['C']
    return values


def stats(path):
    path = Path(path)
    lengths, genes, transcripts, proteins = [], set(), set(), {}
    sequence, name = [], None
    for line in lines(path / 'braker.gtf'):
        columns = line.rstrip().split('\t')
        if len(columns) != 9:
            continue
        if columns[2] == 'gene':
            match = re.search(r'gene_id "([^"]+)"', columns[8])
            identifier = match.group(1) if match else columns[8].strip()
            assert identifier not in genes
            genes.add(identifier)
            length = int(columns[4]) - int(columns[3]) + 1
            assert length > 0
            lengths.append(length)
        if columns[2] == 'transcript':
            match = re.search(r'transcript_id "([^"]+)"', columns[8])
            transcripts.add(match.group(1) if match else columns[8].strip())

    def finish():
        if name is not None:
            translated = ''.join(sequence).rstrip('*')
            assert translated and '*' not in translated, 'Empty protein/internal stop: ' + str(name)
            assert name not in proteins, 'Duplicate protein ID'
            proteins[name] = len(translated)

    for line in lines(path / 'braker.aa'):
        if line.startswith('>'):
            finish()
            name, sequence = line[1:].split()[0], []
        else:
            sequence.append(line.strip())
    finish()
    assert proteins and genes
    assert set(proteins) == transcripts, 'GTF/protein transcript IDs disagree'
    counts, introns, supported_introns, support_ids = Counter(), 0, 0, set()
    rows = csv.DictReader((line for line in lines(path / 'gene_support.tsv') if not line.startswith('#')), delimiter='\t')
    for row in rows:
        identifier = row['transcript_id']
        assert identifier not in support_ids
        support_ids.add(identifier)
        total, supported = int(row['num_introns']), int(row['introns_sup_any'])
        assert 0 <= supported <= total
        category = ('single_exon' if not total else 'all_introns_supported' if total == supported
                    else 'some_introns_supported' if supported else 'no_intron_support')
        counts[category] += 1
        introns += total
        supported_introns += supported
    assert support_ids == set(proteins)
    return dict(genes=len(genes), proteins=len(proteins), mean_gene_bp=statistics.mean(lengths),
                median_gene_bp=statistics.median(lengths), mean_protein_aa=statistics.mean(proteins.values()),
                single_exon_transcripts=counts['single_exon'],
                single_exon_percent=counts['single_exon'] * 100 / len(proteins),
                transcript_support=dict(counts),
                intron_support_percent=supported_introns * 100 / introns if introns else None,
                internal_stops=0, protein_sha256=engine.sha(path / 'braker.aa'))


def choose(rows, target):
    """Preserve the production operational rule, not an accuracy claim."""
    assert rows
    assert all(row['assessments']['arachnida']['N'] == target['N'] for row in rows)
    eligible = [row for row in rows if row['assessments']['arachnida']['C'] >= target['C']]

    def compact(row):
        return (row['annotation']['genes'], row['annotation']['proteins'],
                -row['parameters']['intron_support'], row['label'])

    if eligible:
        return min(eligible, key=compact), ('Genome C reached: minimize gene count, then protein count; '
                'break remaining ties by stricter intron support, then label.')
    best = max(row['assessments']['arachnida']['C'] for row in rows)
    return min([row for row in rows if row['assessments']['arachnida']['C'] == best], key=compact), (
        'No setting reached genome C: maximize C, then minimize model counts using the same tie rules. '
        'Flagged below target; review is required.')


def report(root, rows, target, complete=False):
    directory = root / '03_reports'
    payload = dict(pipeline_version=VERSION, updated=stamp(), all_tests_complete=complete,
                   genome_target=target, trials=rows)
    text = ['# ' + root.name + ' gene-filtering report', '',
            'Status: ' + ('all settings completed' if complete else 'running; no final selection published'), '',
            '## Fixed conditions', '',
            'BRAKER4 ' + engine.MODE.upper() + '; identical upstream predictions/hints; no BUSCO-directed recovery.',
            ('EP: A-H combine AUGUSTUS and raw GeneMark-EP, retaining the high-confidence GeneMark set, '
             'using second-pass protein hints. The native EP baseline uses a different keep-style integration '
             'and is not a strict single-variable control.' if engine.MODE == 'ep' else
             'ETP: AUGUSTUS and GeneMark-ETP candidates; high-confidence training models are retained.'),
            'P=1, E=20, C=1, M=1; e_1=0.1, e_2=0.5, e_3=0.05, e_4=0.2. BRAKER postprocessing is unchanged.',
            'All protein isoforms are assessed independently with the configured compleasm executable: '
            'arachnida_odb12.2 is primary; arthropoda_odb12.2 is auxiliary.',
            'stasto_support=2 disables passing on start/stop support alone; =1 allows this route when both are supported.',
            'Single-exon filtering does not delete every single-exon gene; forced high-confidence retention and other rules still apply.', '',
            'Genome reference C: {:.4f}% ({}/{}).'.format(target['C_percent'], target['C'], target['N']), '',
            '## Settings and results', '',
            '| Group | Intron threshold | Start/stop threshold | Single-exon filter | Genes | Proteins | Mean gene bp | Arachnida C% | S/D/F/M counts | Arthropoda C% | Single-exon transcripts % |',
            '|---|---:|---:|---|---:|---:|---:|---:|---|---:|---:|']
    fields = ['label', 'intron_support', 'stasto_support', 'filter_single_exon_genes', 'genes',
              'proteins', 'mean_gene_bp', 'arachnida_C_percent', 'S', 'D', 'F', 'M',
              'arthropoda_C_percent', 'single_exon_percent', 'intron_support_percent']
    with (directory / 'filtering_results.tsv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(fields)
        for row in rows:
            p, a = row['parameters'], row['annotation']
            primary, auxiliary = row['assessments']['arachnida'], row['assessments']['arthropoda']
            writer.writerow([row['label'], p['intron_support'], p['stasto_support'], p['filter_single_exon_genes'],
                a['genes'], a['proteins'], a['mean_gene_bp'], primary['C_percent'], primary['S'], primary['D'],
                primary['F'], primary['M'], auxiliary['C_percent'], a['single_exon_percent'], a['intron_support_percent']])
            text.append('| {} | {} | {} | {} | {} | {} | {:.1f} | {:.4f} | {}/{}/{}/{} | {:.4f} | {:.2f} |'.format(
                row['label'], p['intron_support'], p['stasto_support'], 'on' if p['filter_single_exon_genes'] else 'off',
                a['genes'], a['proteins'], a['mean_gene_bp'], primary['C_percent'], primary['S'], primary['D'],
                primary['F'], primary['M'], auxiliary['C_percent'], a['single_exon_percent']))
    for row in rows:
        a = row['annotation']
        supported = 'n/a' if a['intron_support_percent'] is None else '{:.3f}%'.format(a['intron_support_percent'])
        text += ['', '### ' + row['label'] + ' additional metrics', '',
                 'Median gene span {:.1f} bp; mean protein length {:.1f} aa; supported introns {}.'.format(
                     a['median_gene_bp'], a['mean_protein_aa'], supported),
                 'Transcript evidence: ' + json.dumps(a['transcript_support']) + '.',
                 'Results: `' + row['result_path'] + '`.']
    if complete:
        winner, reason = choose(rows, target)
        selection = dict(label=winner['label'], reason=reason, parameters=winner['parameters'],
                         C_percent=winner['assessments']['arachnida']['C_percent'],
                         target_met=winner['assessments']['arachnida']['C'] >= target['C'],
                         candidate_labels=[row['label'] for row in rows],
                         limitations='Operational choice within this panel, not a biologically optimal threshold or proof of accuracy.')
        write_json(directory / 'selection.json', selection)
        payload['selection'] = selection
        text += ['', '## Selected setting: ' + winner['label'], '', reason, '',
                 'Completeness takes priority, followed by model-count restraint; this does not maximize accuracy. '
                 'Protein C above genome C does not prove a more accurate annotation.',
                 'Review large increases in single-exon models. All candidates are retained in their trial archives; '
                 'genes are not restored by BUSCO identity.']
    text += ['', '## Interpretation limits', '',
             'Gene length is GTF end-start+1, including introns. D counts duplicated complete markers; '
             'all-isoform inputs affect D, so it is not a direct measure of biological duplication.',
             'No plots or fitted curves are produced. Genome/protein C is compared within the same lineage. '
             'Final output is published only after every setting and its structural/assessment checks succeed.',
             'Version: ' + VERSION + '; updated: ' + stamp()]
    (directory / 'filtering_report.md').write_text('\n'.join(text) + '\n', encoding='utf-8')
    write_json(directory / 'filtering_results.json', payload)
    return choose(rows, target)[0] if complete else None


def collect_case(case, label, parameters):
    result = case / 'results'
    return dict(label=label, parameters=dict(intron_support=float(parameters[0]),
                stasto_support=float(parameters[1]), filter_single_exon_genes=parameters[2]),
                annotation=stats(result / 'annotation'),
                assessments={lineage: assessment(result / ('compleasm_' + lineage), lineage)
                             for lineage in ('arachnida', 'arthropoda')}, result_path=str(result))


def verify_completed(root):
    """Validate all archived files and the published selection, without new assessment."""
    root = engine.no_symlinks(Path(root))
    completion = read_json(root / 'pipeline.complete.json')
    assert completion['version'] == VERSION, 'Completed result is from another pipeline version'
    final = root / '04_final/selected'
    assert completion['final'] == str(final)
    verify_manifest(final)
    for name in REQUIRED_ANNOTATION:
        assert (final / name).is_file() and (final / name).stat().st_size > 0, 'Missing final annotation: ' + name
    provenance = read_json(root / '01_braker/provenance.json')
    archive = root / '01_braker/upstream'
    assert provenance['upstream_sha256'], 'Missing upstream fingerprints'
    for relative, expected in provenance['upstream_sha256'].items():
        rel = Path(relative)
        assert not rel.is_absolute() and '..' not in rel.parts
        assert engine.sha(engine.no_symlinks(archive / rel)) == expected
    assert engine.sha(root / '01_braker/results/annotation/braker.aa') == provenance['baseline_sha256']
    report_data = read_json(root / '03_reports/filtering_results.json')
    assert report_data['all_tests_complete'] is True
    rows = report_data['trials']
    assert [row['label'] for row in rows] == [row[0] for row in PANEL], 'Incomplete/changed threshold panel'
    for row, (_, intron, stasto, single) in zip(rows, PANEL):
        case_result = root / '02_filtering' / row['label'] / 'results'
        verify_manifest(case_result)
        assert row['result_path'] == str(case_result)
        if row['label'] != 'baseline':
            expected = dict(intron_support=float(intron), stasto_support=float(stasto), filter_single_exon_genes=single)
            assert row['parameters'] == expected, 'Threshold panel parameters changed'
            case = case_result.parent
            archived = read_json(case / 'parameters.json')
            assert all(archived[key] == value for key, value in expected.items())
            assert read_json(case / 'upstream_sha256.json') == provenance['upstream_sha256']
            assert read_json(case / 'complete.json')['results'] == str(case_result)
        else:
            assert row['parameters']['intron_support'] == 1.0 and row['parameters']['stasto_support'] == 2
            assert type(row['parameters']['filter_single_exon_genes']) is bool
        assert stats(case_result / 'annotation') == row['annotation'], 'Archived annotation metrics changed'
        assert {lin: assessment(case_result / ('compleasm_' + lin), lin) for lin in ('arachnida', 'arthropoda')} == row['assessments']
    winner, _ = choose(rows, report_data['genome_target'])
    selection = read_json(root / '03_reports/selection.json')
    assert completion['selected'] == winner['label'] == selection['label']
    assert selection == read_json(final / 'selection.json') == report_data['selection']
    assert read_json(final / 'parameters.json') == winner['parameters']
    for name in REQUIRED_ANNOTATION:
        assert engine.sha(final / name) == engine.sha(Path(winner['result_path']) / 'annotation' / name), 'Final file differs from selected trial'
    assert (root / '03_reports/filtering_report.md').is_file()
    return completion


def validate_roots(cfg):
    sample = cfg['sample']
    assert re.fullmatch(r'[A-Za-z0-9_-]+', sample) and sample not in ('.', '..')
    root = engine.no_symlinks(Path(cfg['output_root']))
    temporary = engine.no_symlinks(Path(cfg['temporary_root']))
    assert root not in {Path('/'), Path.home().resolve()} and temporary not in {Path('/'), Path.home().resolve()}
    assert temporary.name == sample
    assert root != temporary and root not in temporary.parents and temporary not in root.parents, 'Result and temporary roots must be disjoint'
    return root, temporary


def own_temporary(temporary, root):
    owner = temporary / '.pipeline-owned'
    engine.no_symlinks(owner)
    if owner.exists():
        assert owner.read_text().strip() == str(root), 'Temporary pipeline ownership mismatch'
    else:
        if temporary.exists() and any(temporary.iterdir()):
            assert {p.name for p in temporary.iterdir()} <= {'.gene-driver-owned', 'qc', 'container_tmp'}, 'Unowned nonempty temporary directory'
            driver = engine.no_symlinks(temporary / '.gene-driver-owned')
            assert driver.is_file() and driver.read_text().strip() == str(root), 'Driver ownership mismatch'
            for name in ('qc', 'container_tmp'):
                p = engine.no_symlinks(temporary / name)
                assert not p.exists() or p.is_dir()
        temporary.mkdir(parents=True, exist_ok=True)
        owner.write_text(str(root) + '\n')
    return owner


def run(cfg):
    if sys.flags.optimize:
        raise RuntimeError('Do not use Python optimization: runtime assertions must remain enabled')
    root, temporary = validate_roots(cfg)
    sample, mode = cfg['sample'], cfg.get('mode', 'etp')
    engine.set_mode(mode)
    threads = int(cfg.get('threads', 56))
    assert 1 <= threads <= 56
    template = configparser.ConfigParser()
    assert template.read(cfg['braker_template']), 'BRAKER template missing'
    assert template.getint('PARAMS', 'run_best_by_compleasm') == 0
    assert template.getint('PARAMS', 'use_compleasm_hints') == 1, 'Preserve native mandatory compleasm hints'
    assert template.getint('PARAMS', 'skip_busco') == 1, 'Keep the agreed BUSCO-evaluation setting'
    for name in ('00_config', '01_braker', '02_filtering', '03_reports', '04_final', 'logs'):
        directory = engine.no_symlinks(root / name)
        directory.mkdir(parents=True, exist_ok=True)
    configdir = root / '00_config'
    saved = configdir / 'pipeline.json'
    if saved.exists():
        assert read_json(saved) == cfg, 'Changed config requires a new output root'
    else:
        write_json(saved, cfg)
    if (root / 'pipeline.complete.json').exists():
        verify_completed(root)
        print('Previously completed archive verified; no analysis repeated.', flush=True)
        return
    owner = own_temporary(temporary, root)
    # Validate scratch ownership before a runtime cache directory can create
    # anything below that scratch root. Completed reuse needs no new cache.
    engine.configure_runtime(cfg['runtime'], template['containers']['braker3_image'])
    copy(Path(cfg['samples_csv']), configdir / 'samples.csv')
    samples = list(csv.DictReader(lines(configdir / 'samples.csv')))
    assert len(samples) == 1 and samples[0]['sample_name'] == sample, 'One genome per internal manifest'
    assert samples[0]['protein_fasta'], 'Protein evidence required'
    if mode == 'etp':
        assert samples[0]['fastq_r1'] and samples[0]['fastq_r2'], 'ETP RNA inputs required'
    else:
        assert not any(samples[0].get(key, '') for key in ['fastq_r1', 'fastq_r2', 'bam_files', 'sra_ids', 'isoseq_bam', 'isoseq_fastq']), 'EP must not silently use RNA'
    template['paths']['samples_file'] = str(configdir / 'samples.csv')
    template['paths']['augustus_config_path'] = str(temporary / '01_braker/augustus_config')
    template['paths']['compleasm_download_path'] = engine.RUNTIME['compleasm_db']
    template['SLURM_ARGS']['cpus_per_task'] = str(threads)
    # memory_mb is the aggregate scheduler budget; preserve a smaller template
    # per-rule request, but do not leave requests above the configured budget.
    rule_memory = template.getint('SLURM_ARGS', 'mem_of_node', fallback=engine.RUNTIME['memory_mb'])
    template['SLURM_ARGS']['mem_of_node'] = str(min(rule_memory, engine.RUNTIME['memory_mb']))
    with (configdir / 'config.ini').open('w') as stream:
        template.write(stream)
    workflow = Path(cfg['workflow'])
    assert workflow.is_absolute() and workflow.is_file()
    for key in ('braker3_image', 'agat_image', 'busco_image'):
        image = Path(template['containers'][key])
        assert image.is_absolute() and image.is_file(), 'Missing local container image: ' + key
    for lineage in ('arachnida', 'arthropoda'):
        database = Path(engine.RUNTIME['compleasm_db']) / (lineage + '_odb12.2')
        for name in ('lengths_cutoff', 'scores_cutoff', 'refseq_db.faa.gz'):
            assert (database / name).is_file(), 'Missing compleasm cache file: ' + str(database / name)
        assert (database / 'hmms').is_dir() and any((database / 'hmms').iterdir()), 'Missing HMM cache'
    target = assessment(Path(cfg['genome_compleasm_summary']).parent, 'arachnida')
    copy(Path(cfg['genome_compleasm_summary']), configdir / 'genome_compleasm.summary.txt')
    if cfg.get('reuse_braker_output'):
        base = Path(cfg['reuse_braker_output'])
        assert base.is_absolute() and Path(cfg['reuse_braker_complete']).is_file()
    else:
        work = engine.no_symlinks(temporary / '01_braker')
        work.mkdir(exist_ok=True)
        env = engine.ENV.copy()
        env['BRAKER4_CONFIG'] = str(configdir / 'config.ini')
        command = [engine.RUNTIME['snakemake'], '--snakefile', str(workflow), '--directory', str(work),
                   '--cores', str(threads)] + engine.snakemake_runtime_arguments() + ['--printshellcmds', '--rerun-incomplete']
        engine.call(command + ['--dry-run'], root / 'logs/braker.dryrun.log', env)
        engine.call(command, root / 'logs/braker.log', env)
        base = work / 'output' / sample
        assert (base / 'results/.done').is_file()
    assert (base / 'braker.aa').is_file()
    tsebra_log = (base / 'tsebra.log').read_text()
    assert 'single-exon gene filter active' in tsebra_log or 'single-exon gene filter inactive' in tsebra_log
    baseline_filter = mode == 'etp' and 'single-exon gene filter active' in tsebra_log
    panel = [('baseline', 1.0, 2, baseline_filter)] + PANEL[1:]
    if cfg.get('baseline_protein_sha256'):
        assert engine.sha(base / 'braker.aa') == cfg['baseline_protein_sha256']
    if not cfg.get('reuse_braker_output'):
        for relative in ('logs', 'benchmarks', 'augustus_config'):
            source, destination = temporary / '01_braker' / relative, root / '01_braker' / relative
            if source.is_dir():
                engine.copy_tree(source, destination)
        for name in ('augustus.hints.gtf', 'augustus.hints.aa', 'augustus.hints.fixed.gtf',
                     'tsebra.log', 'software_versions.tsv', 'hintsfile.gff', 'hintsfile_iter2.gff'):
            if (base / name).is_file():
                copy(base / name, root / '01_braker/key_intermediates' / name)
    archive = root / '01_braker/upstream'
    for relative in engine.INPUTS:
        copy(base / relative, archive / relative)
    for name in KEEP:
        if (base / name).is_file():
            copy(base / name, root / '01_braker/results/annotation' / name)
    write_json(root / '01_braker/provenance.json', dict(source=str(base), mode=mode, workflow=str(workflow),
               workflow_sha256=engine.sha(workflow), baseline_sha256=engine.sha(base / 'braker.aa'),
               upstream_sha256={relative: engine.sha(archive / relative) for relative in engine.INPUTS}))
    engine.configure(root / '02_filtering', temporary / '02_filtering', archive, configdir, workflow, sample, threads)
    baseline = root / '02_filtering/baseline'
    baseline.mkdir(exist_ok=True)
    # Materialize baseline archives: result checksum checks never follow links.
    engine.copy_tree(root / '01_braker/results/annotation', baseline / 'results/annotation')
    for lineage in ('arachnida', 'arthropoda'):
        destination = baseline / 'results' / ('compleasm_' + lineage)
        source = cfg.get('baseline_assessments', {}).get(lineage)
        if source:
            engine.copy_tree(Path(source), destination)
        elif not (destination / 'summary.txt').exists():
            engine.call([engine.RUNTIME['compleasm'], 'protein', '-p', str(base / 'braker.aa'), '-l', lineage,
                         '--odb', 'odb12.2', '-t', str(threads), '-L', engine.RUNTIME['compleasm_db'],
                         '-o', str(destination)], root / 'logs' / ('baseline_' + lineage + '.log'))
    engine.write_manifest(baseline / 'results')
    rows = []
    for label, intron, stasto, single in panel:
        case = root / '02_filtering' / label
        if label != 'baseline':
            old = cfg.get('reuse_trials', {}).get(label)
            if old:
                old = Path(old)
                engine.verify_case(old, label, intron, stasto, single)
                case.mkdir(exist_ok=True)
                engine.copy_tree(old / 'results', case / 'results')
                engine.copy_tree(old / 'logs', case / 'logs')
                for name in ('parameters.json', 'tsebra.cfg', 'complete.json', 'upstream_sha256.json'):
                    copy(old / name, case / name)
                completion = json.loads((case / 'complete.json').read_text())
                completion['results'] = str(case / 'results')
                write_json(case / 'complete.json', completion)
                write_json(case / 'reused_from.json', dict(path=str(old), verified_upstream=True))
            else:
                engine.one(label, str(intron), str(stasto), single)
        verify_manifest(case / 'results')
        row = collect_case(case, label, (intron, stasto, single))
        assert row['assessments']['arachnida']['N'] == target['N'], 'Reference marker count mismatch'
        rows.append(row)
        write_json(case / 'metrics.json', row)
        report(root, rows, target)
        print(stamp(), label, 'VALIDATED C={:.4f}% genes={}'.format(
            row['assessments']['arachnida']['C_percent'], row['annotation']['genes']), flush=True)
        work = temporary / '02_filtering' / label
        if not cfg.get('keep_trial_work', False) and work.exists() and label not in cfg.get('reuse_trials', {}):
            assert owner.read_text().strip() == str(root)
            engine.no_symlinks(work)
            assert work.parent == temporary / '02_filtering'
            marker = engine.no_symlinks(work / '.filter-trial-owned')
            assert json.loads(marker.read_text()) == engine.trial_owner(work, label), 'Trial ownership mismatch'
            engine.verify_case(case, label, intron, stasto, single)
            assert shutil.rmtree.avoids_symlink_attacks, 'fd-based safe recursive removal required'
            shutil.rmtree(work)
            write_json(case / 'cleanup.json', dict(removed_only_temporary_workspace=str(work),
                       time=stamp(), results_retained=str(case / 'results')))
            completion = json.loads((case / 'complete.json').read_text())
            completion.update(work_retained=None, temporary_work_removed=str(work))
            write_json(case / 'complete.json', completion)
    assert len(rows) == len(PANEL)
    winner = report(root, rows, target, True)
    stage, final = root / '04_final/selected.partial', root / '04_final/selected'
    assert not stage.exists() and not stage.is_symlink() and not final.exists() and not final.is_symlink(), 'Inspect existing final publication before retry'
    stage.mkdir()
    source = Path(winner['result_path'])
    for name in KEEP:
        if (source / 'annotation' / name).is_file():
            copy(source / 'annotation' / name, stage / name)
    for lineage in ('arachnida', 'arthropoda'):
        for name in ('summary.txt', 'full_table.tsv'):
            copy(source / ('compleasm_' + lineage) / name, stage / ('compleasm_' + lineage) / name)
    for name in REQUIRED_ANNOTATION:
        assert (stage / name).is_file() and (stage / name).stat().st_size > 0
    copy(root / '03_reports/selection.json', stage / 'selection.json')
    parameters = winner['parameters']
    (stage / 'tsebra.cfg').write_text('P 1\nE 20\nC 1\nM 1\nintron_support {}\nstasto_support {}\ne_1 0.1\ne_2 0.5\ne_3 0.05\ne_4 0.2\n'.format(
        parameters['intron_support'], parameters['stasto_support']))
    write_json(stage / 'parameters.json', parameters)
    engine.write_manifest(stage)
    stage.rename(final)
    verify_manifest(final)
    assert engine.sha(base / 'braker.aa') == json.loads((root / '01_braker/provenance.json').read_text())['baseline_sha256']
    write_json(root / 'pipeline.complete.json', dict(time=stamp(), selected=winner['label'], final=str(final), version=VERSION))
    verify_completed(root)
    write_json(root / 'logs/current.json', dict(stage='all_complete', selected=winner['label'], time=stamp()))
    print(stamp(), 'ALL_COMPLETE selected=' + winner['label'], flush=True)


def main():
    import fcntl
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    args = parser.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    root, _ = validate_roots(cfg)
    logs = engine.no_symlinks(root / 'logs')
    logs.mkdir(parents=True, exist_ok=True)
    lock_path = engine.no_symlinks(logs / 'pipeline.lock')
    with lock_path.open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        # Keep this CLI lease alive in scientific descendants as well.
        existing = os.environ.get('GENEANNOT_LEASE_FDS', '')
        os.environ['GENEANNOT_LEASE_FDS'] = ','.join(filter(None, (existing, str(lock.fileno()))))
        try:
            run(cfg)
        except BaseException as exc:
            write_json(root / 'pipeline.failed.json', dict(time=stamp(), error=repr(exc),
                       final_not_published=not (root / 'pipeline.complete.json').exists()))
            raise
        finally:
            if existing:
                os.environ['GENEANNOT_LEASE_FDS'] = existing
            else:
                os.environ.pop('GENEANNOT_LEASE_FDS', None)


if __name__ == '__main__':
    main()
