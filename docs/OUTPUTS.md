# Output contracts and cleanup

Each sample is stored under `OUTPUT_ROOT/results/SAMPLE/`:

| Location | Retained content |
| --- | --- |
| `00_config` | Input/source hashes, effective config/manifest, supplied genome summary, fastp JSON/HTML and QC summary (ETP). |
| `01_braker/upstream` | Matching genome/index; shared AUGUSTUS, GeneMark candidates, hints and mode-specific evidence used by every trial. |
| `01_braker/results/annotation` | Native baseline GTF/GFF3, proteins, CDS and support records. |
| `01_braker/key_intermediates` | Selected prediction/hint/version files when present. |
| `01_braker/{logs,benchmarks,augustus_config}` | Native run records/training configuration when present. |
| `02_filtering/baseline/results` | Baseline annotation and both protein-compleasm outputs, with checksums. |
| `02_filtering/A` through `H` | Parameters, TSEBRA config, upstream hashes, logs, metrics and checked results for each profile. |
| `03_reports` | Text/TSV/JSON results and automatic selection rationale. |
| `04_final/selected` | Published GFF3/GTF, all-isoform proteins, CDS, evidence summary, profile parameters, lineage summaries/tables, checksums. |
| `logs` | Driver, prediction, assessment and progress logs. |
| `pipeline.complete.json` | Written only after all nine profiles pass and a final selection is verified. |

Final filenames remain `braker.gff3`, `braker.gtf`, `braker.aa`,
`braker.codingseq` and `gene_support.tsv`. Transcript/protein IDs are native;
there is no cross-species ID renaming or longest-isoform extraction. Genome
sequence IDs must match the annotation. Native header cleaning can remove
descriptions; the matching assembly is `01_braker/upstream/genome.fa`.

Reports include gene/protein counts, mean/median gene span (introns included),
mean protein length, single-exon and intron-support metrics, primary and
auxiliary completeness, settings and target status. Protein completeness is
calculated on **all isoforms**, not a longest-isoform subset. This affects D.
Filtering reports are text only; the unchanged native backend may generate its
own graphical/HTML files while predicting, which are not the panel's reports.

## Safe success-only cleanup

Per-profile results are copied from scratch to permanent storage. Exact
SHA-256 manifests reject missing, altered, extra and symlinked result files.
Ownership markers and identical shared evidence are checked before removing
the profile scratch. After final validation and unchanged input/source hashes,
the public driver removes the owned sample scratch when `cleanup_work=true`.
When every sample succeeds and the work root contains only its own marker,
that empty work root is removed as well. Persistent output is not deleted.

Large alignments, cleaned FASTQs and temporary tool products are not retained
by default. Raw RNA, reference proteins, original genome, databases, external
images and caches are never cleanup targets. ETP fastp reports and clean-read
hashes are retained. If later read-level manual curation needs BAMs, rebuild
them from the original reads; they are not promised archive outputs here.

Failure leaves diagnostic workspaces in place. The wrapper does not retry a
failed sample into the same partially populated workspace or relabel it as
complete. Use a fresh dedicated batch output/work pair after diagnosis, or
deliberately remove that sample's failed output and matching owned scratch
after ensuring no job is using them. Never erase a broad project directory.

Verified skipping requires identical manifests, paths, runtime configuration,
input hashes and workflow source hashes. Completed archives embed absolute
paths in provenance: copying the entire archive elsewhere preserves files but
does not transparently support in-place resume from the new path. A checksum
manifest detects corruption; it is not a digital signature or proof against
an attacker who edits both content and checksums.
