# Citations and version reporting

Report the software actually used, its version or source identifier, and the
evidence supplied for each genome. This project wraps existing predictors; it
does not replace their scientific attribution. License information is in
[`THIRD_PARTY_NOTICES.md`](../THIRD_PARTY_NOTICES.md).

## Distinguish the workflow from its backend

- **BRAKER4 workflow source:** Henning Krall and Katharina J. Hoff,
  [BRAKER4](https://github.com/Gaius-Augustus/BRAKER4). Cite the repository and
  the exact snapshot identifier recorded by this release. Do not invent a
  BRAKER4 paper, release number or upstream commit when provenance supplies
  only source-file checksums.
- **Backend container:** the tested image reference is
  `teambraker/braker3:v3.1.1`. This describes the underlying tool container;
  it does not make the frontend "BRAKER4 version 3.1.1". Record a local image
  checksum or immutable digest when available, and retain actual per-tool
  version output. Unmeasured component versions must remain unspecified.
- **This wrapper:** report the gene-annotation-v1 release or commit and the
  selected filtering profile. An automatically selected profile is a
  reproducible decision within the tested panel, not a new prediction method.

## Core methods

Use the references appropriate to the evidence mode. EP and ETP have different
evidence conditions; do not describe EP as having RNA evidence.

| Component or mode | Reference |
| --- | --- |
| ETP / BRAKER3 prediction method | Gabriel L., Brůna T., Hoff K. J., Ebel M., Lomsadze A., Borodovsky M. and Stanke M. (2024). **BRAKER3: fully automated genome annotation using RNA-seq and protein evidence with GeneMark-ETP, AUGUSTUS, and TSEBRA.** Genome Research 34:769–777. [doi:10.1101/gr.278090.123](https://doi.org/10.1101/gr.278090.123). |
| EP / BRAKER2 prediction method | Brůna T., Hoff K. J., Lomsadze A., Stanke M. and Borodovsky M. (2021). **BRAKER2: automatic eukaryotic genome annotation with GeneMark-EP+ and AUGUSTUS supported by a protein database.** NAR Genomics and Bioinformatics 3:lqaa108. [doi:10.1093/nargab/lqaa108](https://doi.org/10.1093/nargab/lqaa108). |
| GeneMark-ETP, when used | Brůna T., Lomsadze A. and Borodovsky M. (2024). **GeneMark-ETP significantly improves the accuracy of automatic annotation of large eukaryotic genomes.** Genome Research 34:757–768. [doi:10.1101/gr.278373.123](https://doi.org/10.1101/gr.278373.123). |
| GeneMark-EP+, when used | Brůna T., Lomsadze A. and Borodovsky M. (2020). **GeneMark-EP+: eukaryotic gene prediction with self-training in the space of genes and proteins.** NAR Genomics and Bioinformatics 2:lqaa026. [doi:10.1093/nargab/lqaa026](https://doi.org/10.1093/nargab/lqaa026). |
| AUGUSTUS | Stanke M., Diekhans M., Baertsch R. and Haussler D. (2008). **Using native and syntenically mapped cDNA alignments to improve de novo gene finding.** Bioinformatics 24:637–644. [doi:10.1093/bioinformatics/btn013](https://doi.org/10.1093/bioinformatics/btn013). |
| TSEBRA | Gabriel L., Hoff K. J., Brůna T., Borodovsky M. and Stanke M. (2021). **TSEBRA: transcript selector for BRAKER.** BMC Bioinformatics 22:566. [doi:10.1186/s12859-021-04482-0](https://doi.org/10.1186/s12859-021-04482-0). |
| compleasm | Huang N. and Li H. (2023). **compleasm: a faster and more accurate reimplementation of BUSCO.** Bioinformatics 39:btad595. [doi:10.1093/bioinformatics/btad595](https://doi.org/10.1093/bioinformatics/btad595). The tested profile uses [compleasm 0.2.9](https://github.com/huangnengCSU/compleasm/releases/tag/v0.2.9). |
| fastp, if this wrapper performs read preprocessing | Chen S., Zhou Y., Chen Y. and Gu J. (2018). **fastp: an ultra-fast all-in-one FASTQ preprocessor.** Bioinformatics 34:i884–i890. [doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560). The tested version is [0.23.3](https://github.com/OpenGene/fastp/releases/tag/v0.23.3), not fastp 1.0. |
| AGAT, for annotation conversion | Dainat J. **AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format.** Zenodo. [doi:10.5281/zenodo.3552717](https://doi.org/10.5281/zenodo.3552717). State the installed version; the tested container reference is 1.4.1. See the [upstream citation instructions](https://github.com/NBISweden/AGAT#how-to-cite). |

## Alignment and supporting tools

Also cite tools actually invoked by the backend, including HISAT2, SAMtools,
DIAMOND, StringTie, ProtHint, miniprot and HMMER where applicable. The upstream
workflow's `report_citations.sh` and generated run-specific citation files are
useful starting points; inspect them against the run logs before publication.
Do not automatically cite optional methods that were disabled.

For the common ETP path, relevant supporting references include:

- Kim D., Paggi J. M., Park C., Bennett C. and Salzberg S. L. (2019). HISAT2 and
  HISAT-genotype. [doi:10.1038/s41587-019-0201-4](https://doi.org/10.1038/s41587-019-0201-4).
- Danecek P. et al. (2021). Twelve years of SAMtools and BCFtools.
  [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008).
- Buchfink B., Reuter K. and Drost H.-G. (2021). Sensitive protein alignments at
  tree-of-life scale using DIAMOND.
  [doi:10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x).
- Kovaka S. et al. (2019). Transcriptome assembly from long-read RNA-seq
  alignments with StringTie2.
  [doi:10.1186/s13059-019-1910-1](https://doi.org/10.1186/s13059-019-1910-1).

If a run reuses user-prepared alignments or cleaned reads, identify their
preparation tools separately rather than implying that this wrapper generated
those inputs. Cite reference-protein collections and lineage databases with
their actual release identifiers and source publications.

## Completeness reporting

The tested selection profile evaluates final proteins with
`arachnida_odb12.2` as the primary lineage and `arthropoda_odb12.2` as an
auxiliary lineage. Record the exact dataset and marker count. Compare genome
and protein completeness only under a matched primary dataset.

Native BRAKER compleasm-derived training hints remain enabled, while
BUSCO-directed post-merge recovery is disabled. Therefore completeness is not
an independent estimate of annotation accuracy. A high C value, or C exceeding
the genome's value, does not by itself establish a biologically correct gene
set. Multiple isoforms can also affect duplicated-marker counts.
