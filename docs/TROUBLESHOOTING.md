# Troubleshooting

- **No lowercase bases:** supply real soft-masked FASTA, not a hard-masked
  or fully uppercase assembly. Lowercase header text does not count.
- **Preflight reports missing lineage files:** finish preparing the correct
  odb12.2 caches outside a run. Do not silently use odb10, odb12 or a different
  organism lineage against the recorded target.
- **GeneMark licensing/key error:** inspect the installed backend's terms and
  setup instructions. No key is supplied by this repository. Do not post keys
  in issues or configuration uploads.
- **Very low RNA mapping:** confirm species, mate pairing, RNA library type,
  contamination, read quality and assembly origin. Check native HISAT2 logs.
  This release does not silently swap RNA species or fall back from ETP to EP.
- **Below-target C:** read all nine profile results. Candidate generation may
  be incomplete; threshold relaxation cannot recover genes never predicted.
  Consider better evidence/assembly or separate model curation rather than
  treating additional relaxed predictions as automatically correct.
- **Too many genes/single-exon models:** inspect repeat overlap, model lengths,
  evidence and representative loci. The retained native postprocessing is
  not a new comprehensive TE-overlap filter. Auto selection balances marker
  recovery and model counts; it does not replace biological QC.
- **Unowned/nonempty workspace:** use new dedicated directories. Do not copy
  marker files into arbitrary directories to defeat ownership checks.
- **Incomplete sample exists:** diagnose its logs and preserve useful evidence;
  then retry with fresh roots or explicitly remove only its failed owned
  sample directories. Other samples can continue, up to the failure limit.
- **Completed file changed:** verification intentionally fails. Do not edit
  archived results in place; make downstream copies for curation.
- **CPU/memory error:** allocate enough CPUs and RAM and align config with that
  allocation. Scheduler budgets are not OS memory enforcement. Do not run two
  drivers against the same managed roots.
- **Unsupported path characters:** relocate/link inputs into shell-safe ASCII
  absolute paths. Managed output/work roots cannot be symlinks.

Please report the package commit, versions, evidence mode, failed stage and
redacted log excerpt. Do not upload real read data or entire private configs.
