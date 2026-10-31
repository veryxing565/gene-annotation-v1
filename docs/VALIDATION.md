# Validation scope

This release separates software validation from biological benchmarking.

The `unittest` suite exercises exact profile definitions, selection and tie
rules, lineage-count mismatches, checksum tampering, extra-file rejection,
symlink/ownership guards, manifest parsing, paired RNA requirements, masking
case checks, fingerprints and resource checks. Linux-only integration tests
run the **actual public driver and core through all nine profiles** with
synthetic outputs substituted only at external tool execution boundaries.
They cover EP, multiple ETP read pairs, mixed-mode serial batches, completed
skip verification, modified final files and failure-limit behavior.

These tests do not measure gene-prediction sensitivity, specificity, RNA
mapping accuracy or biological correctness. The source deployment's use on
real genomes is provenance, not a new benchmark of this portable wrapper.
An EP and ETP Snakemake dry-run can additionally validate the bundled rule
graph without launching predictors. Any recorded dry-run success means DAG
construction only, not a completed real annotation.

GitHub Actions runs synthetic/unit tests on Linux Python 3.11 and validates
source archive creation. Containers, real research data and licensed
GeneMark components are not required for those tests. The optional
`tests/dryrun_workflow.py` uses a locally installed Snakemake/pandas to construct
the EP and ETP DAGs with tiny dummy inputs; no biological commands run.

Before using a new environment at scale, run a small real biological dataset
with your licensed backend images and cached databases. Record image/dataset
checksums and inspect evidence support. No release should claim biological
equivalence merely because mocked integration and dry-runs passed.
