# Dependencies and preparation

## Platform and resource planning

Supported runtime: Linux x86-64 with Python 3.11 and a working Apptainer setup.
Windows is suitable for editing/downloading this source, not running the
predictors. WSL2 requires its own working Linux/Apptainer configuration. The
wrapper uses POSIX locks, process groups, CPU affinity and fd-safe cleanup.
Never use Python `-O` or `PYTHONOPTIMIZE`: core integrity checks use assertions.

Allocate CPUs in your scheduler first. `threads` is 1–56 and cannot exceed
the allocation's CPU affinity. Samples and filtering profiles run serially;
native rules may parallelize within that allocation. `memory_mb` is a
Snakemake aggregate scheduling budget, **not an operating-system memory cap**;
use your scheduler/cgroup for an enforced limit. `rule_memory_mb` is the
configured native per-rule request, while some native rules have fixed requests.
Provide sufficient memory for the genome; the example 240000/120000 MB is
not a universal recommendation or a guarantee of fit.

Temporary disk demand depends on genome size and RNA depth. The default
free-space checks (100 GiB work, 20 GiB output) are only starting thresholds,
not estimates or reservations. If both roots share a filesystem they share
capacity. Place scratch on fast local storage and permanent output on the
project disk; ensure the latter is mounted before starting. Inputs and all
dependencies must be outside the dedicated scratch/output roots.

## Tested versions, not automatically upgraded

| Component | Recorded version |
| --- | --- |
| Host Python | 3.11.15 |
| Snakemake | 8.18.2 |
| Apptainer | 1.5.2 |
| fastp | 0.23.3 |
| Independent panel compleasm | 0.2.9 |
| Bundled BRAKER4 frontend header | 0.5.0-beta; exact source hashes recorded |
| Tool-container reference | `teambraker/braker3:v3.1.1` |
| AUGUSTUS inside recorded backend | 3.5.0, commit 220e5a6 |
| GeneMark-ETP | 1.02, commit aacc025 |
| TSEBRA | commit 23ef205 |
| Native compleasm inside backend | 0.2.8 |
| HISAT2 / StringTie / DIAMOND | 2.2.1 / 2.2.1 / 2.0.15 |
| SAMtools / miniprot | 1.13 / 0.12-r237 |
| AGAT image | 1.4.1 |
| BUSCO image, optional evaluation disabled | 6.1.0 |

These versions were recorded from the source deployment; they are not a
claim that every current upstream release is compatible. No upstream git
commit was available for the full frontend snapshot, and the original SIF
digest was not measured. Image tags alone are not immutable identifiers.

Install host packages in an isolated environment. `requirements.txt` is for
Python packages only; `environment.yml` is an optional recipe, not a solved
lockfile. Apptainer may need administrator configuration on your platform.
See official [Snakemake](https://snakemake.readthedocs.io/),
[Apptainer](https://apptainer.org/docs/user/latest/),
[fastp](https://github.com/OpenGene/fastp) and
[compleasm](https://github.com/huangnengCSU/compleasm) installation instructions.

## Local images

Obtain lawfully prepared compatible SIF images and put their absolute paths
in `containers`. The tested upstream template references the following images;
you may stage them on a network-enabled machine and transfer them to a server:

```bash
apptainer pull braker3-v3.1.1.sif docker://teambraker/braker3:v3.1.1
# Obtain an AGAT 1.4.1 SIF and a BUSCO 6.1.0 SIF from their official distributions.
sha256sum *.sif > images.SHA256SUMS
```

Consult the [BRAKER4](https://github.com/Gaius-Augustus/BRAKER4),
[AGAT](https://github.com/NBISweden/AGAT) and
[BUSCO](https://busco.ezlab.org/) instructions. Do not guess image tags or
assume an image filename proves its contents. Both EP and ETP also require
GeneMark/ProtHint and AUGUSTUS resources as provided by the compatible backend.
Handle licensing and any required key locally; never commit a key here.

All three image paths are required by the preserved frontend configuration,
even though separate BUSCO QC is disabled. `check` probes host executable
versions but does not launch/validate every internal image tool. A normal
run performs the native Snakemake dry-run before prediction. Validate a new
image using a small real biological dataset before a large batch.

## Lineage caches: prepare before running

Use compleasm 0.2.9 to download **both** `arachnida_odb12.2` and
`arthropoda_odb12.2` into a shared library directory. Follow its `download -h`
for the exact interface of the installed version (use `--odb odb12.2`).
The expected flat layout under `databases.compleasm` is:

```text
compleasm/
  arachnida_odb12.2/{lengths_cutoff,scores_cutoff,refseq_db.faa.gz,hmms/...}
  arthropoda_odb12.2/{lengths_cutoff,scores_cutoff,refseq_db.faa.gz,hmms/...}
```

Preserve all other files and cache metadata produced by upstream's downloader,
not just the four items used for minimum preflight checks. Provide the BUSCO
cache path separately or point it to the same directory if compatible with
your local setup; standalone BUSCO is disabled. Native mandatory compleasm
genome/hint stages remain enabled and use the configured local cache.

Prepare the **same assembly's existing genome-mode compleasm 0.2.9 assessment**
with the primary lineage. Supply its `summary.txt` in the samples TSV; the
wrapper uses it as a selection target and does not initiate a separate
replacement genome assessment. It checks lineage name and marker counts,
but cannot prove that a summary really came from that assembly: retain your
own input provenance. Native BRAKER's required genome assessment still runs.

The wrapper does not download missing files. Scientific subprocesses get an
unusable HTTP(S) proxy to prevent accidental proxy use; this is **not a
network-isolation guarantee**. A missing cache should be fixed during setup,
not by enabling network downloads inside an annotation run. Pre-stage files
and enforce network policy at the system level if required.
