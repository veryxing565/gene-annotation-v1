# Upstream source, scope and integration notes

## Source identity

The vendored `workflow/` is a text-source snapshot of
[Gaius-Augustus/BRAKER4](https://github.com/Gaius-Augustus/BRAKER4), including
the workflow rules and helper scripts needed by this release. It is not an
independently reimplemented gene predictor and must not be advertised as an
official upstream BRAKER4 release.

Use the accompanying release provenance and file checksums to identify the
exact source snapshot. A historical local snapshot may already differ from
the current upstream branch. Do not manufacture a tag or commit identifier to
describe it. Compare snapshots before replacing the workflow with an upstream
update; a newer version is not automatically equivalent to the tested one.

The backend image reference `teambraker/braker3:v3.1.1` and the BRAKER4
Snakemake source are separate versioned components. The filename of a SIF
does not establish the versions of TSEBRA, AUGUSTUS, GeneMark or other tools
inside it. Actual images and executables are supplied by the user.

## Deliberately restricted public interface

Only these evidence modes are supported:

- **EP:** protein evidence without RNA evidence.
- **ETP:** protein evidence together with user-provided paired-end RNA-Seq
  reads. The backend constructs its RNA alignments; the public profile does
  not promise a BAM-only or long-read input route.

A pre-existing soft-masked genome is required. Genome sequence IDs and base
case are input data, not labels to standardize during gene annotation. Preserve
the soft mask and make evidence sequence IDs agree with the supplied genome.

This package does not expose internal repeat masking, SRA/VARUS downloads,
IsoSeq/dual/ES/ET-only annotation, ncRNA prediction, OMArk or FANTASIA. Related
Unused helper scripts are excluded; optional rule definitions are retained
for source provenance and shared parsing but are disabled. Do not enable
these features by editing a low-level configuration field; they are outside
the supported interface and may require files not distributed here.

Some rules for other prediction stages must nevertheless remain present:
the upstream Snakefile includes shared ET, ES or protein-hint rules based on
available evidence, even when a particular rule is not selected by an EP/ETP
DAG. Source pruning must preserve both EP and ETP parse/dry-run checks.

## Scientific behavior retained

### Native prediction and postprocessing

The wrapper runs the established BRAKER workflow before applying the filter
panel. Training, candidate generation and native postprocessing are not
replaced by the panel. The original AUGUSTUS/GeneMark predictions, hints,
training configuration, software records and filtering results are retained
as documented by the run's output inventory.

The configuration deliberately keeps:

- `run_best_by_compleasm = 0`: **no BUSCO-directed recovery** of models lost
  during integration.
- `use_compleasm_hints = 1`: the upstream native compleasm-derived hint
  behavior remains enabled; disabling rescue does not disable these hints.
- `skip_busco = 1`: the separate optional BUSCO quality-control rule is not
  used; the wrapper evaluates final proteins with compleasm instead.

The `best_by_compleasm` Snakemake rule remains part of the DAG so the expected
intermediate filenames are preserved. With recovery disabled, the rule copies
the incoming TSEBRA GTF through unchanged and exits before invoking
`best_by_compleasm.py`. That inherited script may remain in the source bundle
for provenance, under its Artistic-2.0/Apache-2.0 notices; its presence is not
permission to enable recovery through the public interface.

### Nine filtering profiles

The baseline is the native workflow result. A–H use the same upstream
candidates and evidence within the same mode; they change only the listed
TSEBRA selection settings before the same postprocessing:

| Profile | `intron_support` | `stasto_support` | Single-exon filter |
| --- | ---: | ---: | --- |
| baseline | 1.00 | 2 | Native mode-dependent behavior |
| A | 0.90 | 2 | enabled |
| B | 0.80 | 2 | enabled |
| C | 0.80 | 1 | enabled |
| D | 0.80 | 1 | disabled |
| E | 0.75 | 2 | enabled |
| F | 0.70 | 2 | enabled |
| G | 0.60 | 2 | enabled |
| H | 0.50 | 2 | enabled |

Other tested TSEBRA settings remain `P=1`, `E=20`, `C=1`, `M=1`,
`e_1=0.1`, `e_2=0.5`, `e_3=0.05`, `e_4=0.2`.

`stasto_support=2` disables acceptance through the start/stop-only support
route; `1` permits the fully supported start/stop condition. Enabling the
single-exon filter does not mean deleting all single-exon genes: high-confidence
keep sets and the other TSEBRA rules still apply.

In ETP, the panel uses AUGUSTUS and GeneMark-ETP candidates, with the training
high-confidence set retained. In EP, it uses AUGUSTUS and GeneMark-EP
candidates, the GeneMark high-confidence set, and second-round protein hints.
The EP raw and good sets receive separate ID namespaces to prevent unrelated
models with colliding IDs from being forced into the retained set.

The EP baseline uses upstream's retain-style integration and is therefore not
a strict threshold-only control for A–H. Comparisons among A–H hold their
candidate set fixed. Do not attribute every EP baseline-versus-panel change
to a single numeric threshold.

### Evaluation and automatic recommendation

Each profile's final proteins are assessed independently by compleasm. The
tested evaluator is 0.2.9, with `arachnida_odb12.2` as the primary dataset and
`arthropoda_odb12.2` as an auxiliary dataset. Databases are prepared by the user;
the run must not silently substitute another release or download a different
dataset. The supplied genome-completeness target must use the same primary
lineage and marker count as the protein evaluations.

Selection is an operational rule over the tested profiles:

1. Among profiles reaching the genome target C, choose the fewest genes.
2. Break ties by the fewest protein sequences, then stricter intron support,
   followed by the deterministic profile label.
3. If no profile reaches the target, choose the highest C with the same
   tie-breakers and flag the target as unmet.

The report records each profile's conditions, gene/protein counts and
completeness, the recommendation and its limitations. The wrapper does not
fit a threshold curve or generate an additional graphical panel report.
Upstream's own backend reporting is a separate component.

This rule does not prove that the selected profile has the best biological
precision or that thresholds transfer optimally to another genome. Native
marker-derived hints and repeated use of marker completeness limit its
independence; use independent evidence when validating annotation accuracy.

## Operational adaptation, not a new biological method

Public configuration and orchestration must use explicit user-provided paths,
bounded CPU resources, separate per-sample temporary workspaces and permanent
result locations. Machine-specific paths, study-specific sample lists,
credentials and research inputs do not belong in the public source tree.

Completion, export verification, resuming, skipping and cleanup are storage
and scheduling behavior. They must not silently change prediction parameters
or treat a failed annotation as complete. Preserve diagnostic workspaces on
failure and remove only owned temporary data after verified publication.

### Public source modifications (2026-10-31)

The site-specific core orchestrator and filter engine were adapted to accept
explicit runtime paths, local containers/databases, configurable resources,
English reporting, ownership/checksum checks and process-group cleanup. The
panel and scientific selection/postprocessing behavior are retained.

Vendored workflow text is retained from the deployment snapshot, with LF
line endings and one portability correction in
`rules/postprocessing/postprocess_augustus.smk`: hint-support analysis uses
the absolute `{script_dir}/analyze_hint_support.py` instead of assuming a
`scripts/` directory in the scratch working directory. This changes script
resolution, not its analysis parameters. Source provenance lists original
and packaged hashes; it does not claim byte-identical copies or an unknown
upstream git commit.

No upstream tools, containers, restricted GeneMark components or license keys
are included. See [`THIRD_PARTY_NOTICES.md`](../THIRD_PARTY_NOTICES.md) for
source-level attribution and [`CITATIONS.md`](CITATIONS.md) for method citations.

Future changes to vendored files should retain the original headers and carry
a dated modification note, a clear description of the change, and updated
source checksums. Do not describe a local compatibility patch as an upstream
bug fix unless it has actually been accepted upstream.
