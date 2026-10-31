# gene-annotation-v1

[![Tests](https://github.com/veryxing565/gene-annotation-v1/actions/workflows/tests.yml/badge.svg)](https://github.com/veryxing565/gene-annotation-v1/actions/workflows/tests.yml)

A portable, **arachnid-focused protein-coding gene annotation workflow**:
paired-end RNA quality control, a preserved BRAKER4 Snakemake source snapshot,
nine TSEBRA selection profiles, independent protein-mode compleasm assessments,
and a reproducible final-profile recommendation. EP (proteins only) and ETP
(RNA-seq plus proteins) are supported.

This is an orchestration package, **not a new predictor or an official BRAKER
release**. Version `1.0.0rc1` is a portable release candidate derived from an
in-use annotation workflow. The public wrapper has software tests and DAG
checks; it has not been benchmarked as a new biological method. See
[validation scope](docs/VALIDATION.md).

## What it does

1. Validate supplied soft-masked genome, reference proteins, paired RNA
   manifests (ETP), local dependencies and an existing genome-completeness target.
2. In ETP, clean each paired RNA run with fastp; in EP, do not fabricate RNA evidence.
3. Run the bundled BRAKER4 workflow with the established settings: native
   compleasm hints enabled, **BUSCO-directed post-merge recovery disabled**.
4. Reuse the same upstream predictions and hints for baseline plus eight
   alternative TSEBRA profiles; apply the same native postprocessing.
5. Assess every profile's proteins with compleasm 0.2.9, using
   `arachnida_odb12.2` (primary) and `arthropoda_odb12.2` (auxiliary).
6. Write English Markdown, TSV and JSON reports, select a profile, verify
   exported results, then remove only owned temporary data on success.

The wrapper runs genomes **serially**, caps CPU affinity at the configured
1–56 CPUs, skips only verified completed results, and stops after three
consecutive sample failures by default. Failed temporary data are preserved.
There are no automatic research-data or database downloads.

## Install dependencies

Use Linux x86-64 with Python 3.11, Snakemake **8.18.2**, fastp **0.23.3**,
compleasm **0.2.9**, Apptainer (tested **1.5.2**) and pandas. Prepare local
BRAKER, AGAT and BUSCO SIF images and both lineage caches separately.

```bash
git clone https://github.com/veryxing565/gene-annotation-v1.git
cd gene-annotation-v1
# In an existing isolated environment; this installs host Python packages only:
python -m pip install -r requirements.txt
# Or use environment.yml as an optional, unsolved conda recipe.
python run_gene_annotation.py --version
```

The underlying tested tool image is `teambraker/braker3:v3.1.1`. That tag is
**not** a BRAKER4 version: the bundled BRAKER4 frontend identifies itself as
`0.5.0-beta`. Native container compleasm is 0.2.8; the separate panel evaluator
is 0.2.9. Read [dependency setup and versions](docs/DEPENDENCIES.md) before use.
Containers, databases, GeneMark software/keys, and research inputs are not
distributed here. Their licenses remain separate.

## Quick start

Copy the configuration and manifests, then replace every placeholder with your
own absolute Linux paths. Paths must use ASCII letters, digits, `_`, `-`, `.`,
`+` and `/` only: no spaces, colons, commas or shell metacharacters.

```bash
cp config/runtime.example.json runtime.local.json
cp examples/samples.tsv samples.local.tsv
cp examples/example_etp.rna.tsv example_etp.rna.local.tsv
# Edit the three files. Do not run placeholder paths unchanged.
python run_gene_annotation.py check --config runtime.local.json --samples samples.local.tsv
python run_gene_annotation.py run --config runtime.local.json --samples samples.local.tsv
```

Use `tmux`, a batch scheduler allocation, or your usual background launcher for
long jobs. Run the same command again to verify and skip completed samples;
incomplete samples are not silently resumed or overwritten. See
[inputs and configuration](docs/USAGE.md) and [failure handling](docs/TROUBLESHOOTING.md).

The genome must be an **uncompressed, genuinely soft-masked FASTA**. Lowercase
bases are preserved, not converted to uppercase. Reference proteins must be an
uncompressed FASTA with meaningful unique IDs. The prior genome compleasm
`summary.txt` must match this assembly and `arachnida_odb12.2`. RNA input is
paired FASTQ/FASTQ.gz, listed as separate runs; same-species evidence is strongly
preferred. Read provenance/species identity are the user's responsibility.

## Main outputs

```text
OUTPUT_ROOT/
  state/batch.json
  results/SAMPLE/
    00_config/                 configuration, input hashes, fastp reports
    01_braker/                 baseline, shared predictions/hints, training logs
    02_filtering/{baseline,A,...,H}/  all profile results and assessments
    03_reports/                filtering_report.md, TSV/JSON, selection.json
    04_final/selected/         chosen GFF3/GTF, proteins, CDS, compleasm summaries
    logs/
    pipeline.complete.json
```

The final annotation is `04_final/selected/braker.gff3`; its proteins are
`braker.aa` and CDS are `braker.codingseq`. Keep the matching
`01_braker/upstream/genome.fa`: upstream can clean FASTA header descriptions.
All protein isoforms are retained. Longest-isoform selection, comparative
genomics ID renaming, GALBA/manual integration and functional annotation are
**not** part of this package. [Output and cleanup details](docs/OUTPUTS.md).

## Selecting a profile is not proving accuracy

The operational rule first considers profiles reaching the prior genome's
complete-marker count (C = S + D), then prefers fewer genes, fewer proteins,
stricter intron support and finally the profile label. If none reaches the
target, it selects the highest C with the same tie-breakers and flags the
target as unmet. More complete markers do not establish gene-model accuracy;
native marker hints also limit independence. EP and ETP are not equivalent
evidence conditions, and all isoforms can inflate duplicated-marker counts.

See the [nine profiles and method boundaries](docs/UPSTREAM_PATCHES.md),
[citations](docs/CITATIONS.md), [source provenance](provenance/source_snapshot.json)
and [third-party notices](THIRD_PARTY_NOTICES.md).

## Development and source archives

```bash
python -m unittest discover -s tests -v
python scripts/build_release.py --output ../gene-v1-source-release
```

Tests use synthetic data and mocked external executables; they do not require
genomes or licensed predictors. Source archives include file checksums.
Original wrapper code is MIT; vendored files retain their own notices. Please
report issues with redacted configuration and versions, not private data.
