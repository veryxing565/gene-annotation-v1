# Input manifests and execution

Run `python run_gene_annotation.py --help`. The public interface has `check`
and `run`; do not invoke the lower-level core directly for normal batches.

## Input files

`samples.tsv` must use Unix LF endings, tabs and exactly these columns:

```text
sample_id  mode  genome  protein_fasta  genome_compleasm_summary  rnaseq_manifest
```

The spaces above are for display; use the actual tab-separated examples.
One row per genome. `sample_id` starts with a letter/digit and contains only
letters/digits, `_` or `-` (maximum 120 characters). Modes are lowercase
`ep` or `etp`. All input paths are absolute. Genome and proteins are
uncompressed FASTA; reads may be FASTQ.gz or FASTQ.

For EP, set `rnaseq_manifest` to `-` or empty. For ETP, give a TSV with exactly
`run_id`, `fastq_r1`, `fastq_r2`. Each row is one paired run, IDs are unique,
mates are distinct nonempty files. Multiple runs remain distinct through QC,
then are supplied to the native frontend as paired lists. Single-end reads,
long-read RNA and BAM-only input are not exposed by this release. Do not mix
species or guessed mates: existence checks do not establish biological identity.

Prefer taxonomically relevant, clean reference proteins. Duplicate sequence
identifiers, genomic contamination, mismatched RNA or shallow/tissue-limited
RNA can harm prediction; high marker completeness does not remove these risks.
The package does not pick reference species or retrieve RNA automatically.

The soft-masked input is checked for lowercase `acgt` in sequence, not headers.
No case conversion is performed. Native preprocessing may strip FASTA header
descriptions and excludes contigs shorter than 10000 bp. Keep the archived
matching genome for downstream coordinate work. Check unusual/duplicate
sequence IDs before submitting a large job.

## Configuration

Edit a copy of `config/runtime.example.json`. Output and work roots must be
dedicated, disjoint, non-symlink directories, separate from this source tree.
They must initially be absent or empty. Reuse requires matching ownership
markers; do not point them at an existing project root or a home directory.
Input symlinks are allowed if their resolved targets are valid and outside
both managed roots. Tool/database/cache paths must also be outside scratch.

| Key | Meaning |
| --- | --- |
| `output_root` | Persistent archive root; never automatically deleted. |
| `work_root` | Local temporary root; owned successful work is removed. |
| `threads` | CPU-affinity and scheduling limit, 1–56. |
| `memory_mb`, `rule_memory_mb` | Aggregate scheduling budget and native per-rule configured request. Not an OS RAM limit. |
| `fastp_threads` | Per-run QC threads, default min(8, threads). |
| `cleanup_work` | Remove completed sample scratch after verification; default true. Individual successful filter scratch is also removed by the core. |
| `max_consecutive_failures` | Default 3; successful samples reset the streak, verified skips do not. |
| `min_free_work_gb`, `min_free_output_gb` | Nonnegative preflight free-space thresholds, measured in GiB. |
| `tools` | Local executable paths, or resolvable executable names. fastp is checked even for an EP-only batch. |
| `containers` | Local SIF paths; backend BRAKER3 tool image, AGAT and BUSCO. |
| `databases` | Local compleasm and BUSCO library roots. |
| `container_cache`, `apptainer_cache` | Persistent caches outside scratch; not removed. |
| `java_home` | Optional Java installation root; otherwise environment default. |

`check` reads inputs and probes versions without creating managed work/output
directories or running predictors. A run hashes all inputs in full before
annotation and again before cleanup; very large read sets add I/O time.
Source hashes, configuration and version-probe output are also retained.

## Fixed biological settings

`config/braker.parameters.ini` records the production profile. Key settings:
`min_contig=10000`, AUGUSTUS optimization on, single-exon training downsampling
on (`lambda=2`), `use_compleasm_hints=1`, `skip_busco=1`, `run_best_by_compleasm=0`,
translation table 1, GC donor 0.001, hinted splice sites `gcag,atac`, AUGUSTUS
chunks 3 Mb with 0.5 Mb overlap. Supplied soft masking is used; internal
masking, ncRNA, OMArk, FANTASIA and minisplice are off.

ETP fastp uses adapter detection for paired reads, poly-G trimming, Q20,
maximum 40% low-quality bases, at most 5 N bases and minimum length 50 bp.
No duplicate removal is requested. The native ETP frontend uses HISAT2 and
StringTie; **STAR is not substituted here**. This package does not add a
new mapping-rate rejection threshold to the production method. Review native
alignment logs and biological evidence before accepting a weak result.

No BUSCO-directed recovery means exactly that post-merge recovery is off;
it does not turn off mandatory native compleasm/hints, high-confidence keep
sets or the native non-BUSCO both-predictors-agree rescue. The nine settings,
EP namespace correction and selection rule are detailed in
[`UPSTREAM_PATCHES.md`](UPSTREAM_PATCHES.md).

## Background jobs and return codes

Run in an allocated scheduler job or terminal multiplexer. There is no daemon,
automatic login, remote transfer or queue scheduler in this repository.
The destination may be on a project disk: archives are written/copied from
local scratch and verified before scratch is removed. No hardcoded server
name, user account or project mount is needed.

Return codes: `0` all selected samples completed/verified skipped, `1` batch
finished with sample failures, `2` configuration/preflight/ownership error,
`3` consecutive-failure limit reached, `128+signal` interruption. Always inspect
`state/batch.json` and the per-sample driver/native logs. Do not delete lock
files to bypass a live process: file existence alone does not mean a lock is held.
