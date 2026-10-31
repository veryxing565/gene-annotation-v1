# Third-party notices

The license at the project root covers this project's original orchestration,
configuration, tests and documentation. It does **not** relicense upstream
software, containers, databases or research inputs. Existing copyright and
license notices in vendored files must remain intact.

## Bundled BRAKER4 workflow source

`workflow/` contains a selected text-source snapshot of the
[BRAKER4 Snakemake workflow](https://github.com/Gaius-Augustus/BRAKER4), developed
by Henning Krall and Katharina J. Hoff. It is not a new gene-prediction engine
and is not represented as an official upstream release.

The BRAKER4 repository's root license is MIT, with copyright (c) 2025 Katharina
Hoff. The original notice is retained in [`workflow/LICENSE`](workflow/LICENSE).
Upstream source identification and file checksums are recorded separately in
the release provenance. Public packaging changes are described in
[`docs/UPSTREAM_PATCHES.md`](docs/UPSTREAM_PATCHES.md).

Some inherited scripts contain more specific notices than the repository-wide
MIT declaration. Those notices are preserved; the project does not claim that
these files have been converted to MIT:

| Bundled source | Attribution and applicable notices |
| --- | --- |
| `workflow/Snakefile` and `workflow/rules/` | BRAKER4 workflow source, covered by its MIT notice unless a more specific notice appears in a file. |
| `workflow/scripts/analyze_hint_support.py`, `both_agree_rescue.py`, `completeness_plot.py`, `gene_set_statistics.py`, `gene_support_summary.py`, `generate_report.py`, `normalize_gtf.py`, `report_citations.sh`, `training_summary.py` | BRAKER4 source retained under its upstream MIT notice; preserve any authorship information in the files. |
| `workflow/scripts/compleasm_to_hints.py` | Katharina J. Hoff, copyright 2023; the file specifies the Artistic License. Preserve its original BRAKER attribution and the Artistic-1.0 text distributed with BRAKER. |
| `workflow/scripts/stringtie2utr.py` | Katharina J. Hoff; the file specifies the Artistic License. Preserve its original BRAKER attribution and the Artistic-1.0 text distributed with BRAKER. |
| `workflow/scripts/getAnnoFastaFromJoingenes.py` | Katharina J. Hoff; derived from the AUGUSTUS script. Preserve its upstream attribution and modification notes. AUGUSTUS scripts are distributed under Artistic-1.0. |
| `workflow/scripts/fix_in_frame_stop_codon_genes.py` | Anica Hoppe and Katharina Hoff; derived from AUGUSTUS. Preserve its attribution and the Artistic-1.0 notice for AUGUSTUS scripts. |
| `workflow/scripts/best_by_compleasm.py` | Katharina J. Hoff, with credit to Neng Huang. The file explicitly identifies Artistic-2.0 and Apache-2.0 portions copied and modified from compleasm. Preserve the file's section-specific notices. The public pipeline disables this script's BUSCO-directed recovery; retaining source is not enabling that feature. |

License texts included with the source:

- [`workflow/LICENSES/Artistic-1.0.txt`](workflow/LICENSES/Artistic-1.0.txt), from
  [BRAKER's LICENSE.TXT](https://github.com/Gaius-Augustus/BRAKER/blob/master/LICENSE.TXT).
- [`workflow/LICENSES/Artistic-2.0.txt`](workflow/LICENSES/Artistic-2.0.txt), from
  [TSEBRA's bin/LICENSE.txt](https://github.com/Gaius-Augustus/TSEBRA/blob/main/bin/LICENSE.txt).
- [`workflow/LICENSES/Apache-2.0.txt`](workflow/LICENSES/Apache-2.0.txt), from
  [compleasm's LICENSE](https://github.com/huangnengCSU/compleasm/blob/main/LICENSE).

These separate texts document file-level provenance; they are not a statement
that every bundled file is offered under every listed license. No upstream
author or institution is claimed to endorse this project.

## Deliberately excluded source and features

This release supports EP and ETP annotation of a **user-supplied, already
soft-masked genome**. It does not provide internal repeat masking, research-data
downloads, VARUS, ncRNA prediction, OMArk or FANTASIA.

In particular, `parseTrfOutput.py` is not bundled: its inherited header names
Georgia Institute of Technology, but this audit did not establish an explicit
file-specific redistribution grant. It is not required when internal masking
is disabled. The internal-masking helper `splitMfasta.pl`, automatic-download
helpers and unused ncRNA/FANTASIA scripts are also excluded. A disabled feature
must not be enabled against a deliberately reduced workflow snapshot.

## External dependencies: supplied by the user

No executable binary, container image, GeneMark key, database, pretrained
model, genome, RNA read file or reference-protein collection is distributed
with this project. Users must obtain and prepare compatible dependencies
themselves, respecting the terms attached to each version they select.

| External component | Official source and license boundary |
| --- | --- |
| BRAKER tool container | The tested backend image reference is `teambraker/braker3:v3.1.1`. This is a **container tag**, not the version of the BRAKER4 Snakemake source. Individual tools inside it retain their own licenses and versions. No image is included here. |
| TSEBRA | [Gaius-Augustus/TSEBRA](https://github.com/Gaius-Augustus/TSEBRA); its Python source specifies [Artistic-2.0](https://github.com/Gaius-Augustus/TSEBRA/blob/main/bin/LICENSE.txt). |
| AUGUSTUS | [Gaius-Augustus/Augustus](https://github.com/Gaius-Augustus/Augustus); source and scripts specify the [Artistic License](https://github.com/Gaius-Augustus/Augustus/blob/master/src/LICENSE.TXT), version 1.0 in the provided text. |
| GeneMark family | [Official licensing/download page](https://genemark.bme.gatech.edu/license_download.cgi) and [GeneMark-ETP](https://github.com/gatech-genemark/GeneMark-ETP). Terms differ by component and distribution: the public ETP pipeline describes CC-BY-NC-SA-4.0, whereas other GeneMark distributions use a restrictive academic agreement. Commercial use may require separate authorization. No key or restricted component may be added to this repository. Some recent distributions no longer require a key; do not assume one licensing model covers every version. |
| compleasm | [huangnengCSU/compleasm](https://github.com/huangnengCSU/compleasm); Apache-2.0, with additional BUSCO attribution in [LICENSE-BUSCO](https://github.com/huangnengCSU/compleasm/blob/main/LICENSE-BUSCO). The tested evaluation profile uses **0.2.9**, not an automatically updated main branch. |
| AGAT | [NBISweden/AGAT](https://github.com/NBISweden/AGAT); GPL-3.0 and its own third-party notices. The tested container reference identifies **1.4.1**; no AGAT code or image is bundled. |
| fastp | [OpenGene/fastp](https://github.com/OpenGene/fastp); MIT. The tested preprocessing profile uses **0.23.3**. No executable is bundled. |
| Other backend tools and Python dependencies | HISAT2, SAMtools, DIAMOND, StringTie, ProtHint, miniprot, HMMER, Snakemake, Apptainer and Python packages retain their upstream terms. Obtain them from their respective projects or a lawfully prepared compatible environment/container. Record actual versions rather than inferring them from a container filename. |

The software license for compleasm or BUSCO does not license their lineage
datasets. Upstream's BUSCO notice identifies the datasets as CC-BY-ND-4.0;
consult the terms supplied with the exact database release. This repository
does not redistribute lineage datasets or reference proteins.

See [`docs/CITATIONS.md`](docs/CITATIONS.md) for scientific attribution. A
citation is not a substitute for complying with a software or data license.
