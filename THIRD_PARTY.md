# Third-party software and data

This repository is MIT-licensed (see [`LICENSE`](LICENSE)). **Nothing
third-party is vendored into it** — every tool is pulled at build or run time by
pinned digest, and every container is built from an upstream base. What follows
is what those tools are, what they are licensed under, and the one case where a
licence needs a word of explanation.

Versions here are what is pinned in `docker/*/Dockerfile` and
`conf/base.config`. `make versions` reports what the built images actually
contain, which is the authority if the two ever disagree.

## Software invoked by the pipeline

| tool | version | licence | how it is used |
|---|---|---|---|
| [Nextflow](https://www.nextflow.io/) | ≥ 23.10.0 | Apache-2.0 | workflow engine; a host requirement, not shipped here |
| [minimap2](https://github.com/lh3/minimap2) | 2.28 | MIT | all read mapping |
| [htslib](https://github.com/samtools/htslib) | 1.21 | MIT | BAM/SAM I/O |
| [samtools](https://github.com/samtools/samtools) | 1.21 | MIT | sorting, indexing, depth |
| [seqkit](https://github.com/shenwei356/seqkit) | 2.8.2 | MIT | FASTA/FASTQ manipulation |
| [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/) | v2 | US Government work, public domain | fetching reference genomes by accession |
| [Debian](https://www.debian.org/) bookworm-slim | 20241111 | mixed, predominantly GPL/LGPL/MIT/BSD | base image for the tools container |
| [Python](https://www.python.org/) | 3.12-slim-bookworm | PSF-2.0 | base image for the analysis container |
| [NumPy](https://numpy.org/) | 2.1.3 | BSD-3-Clause | numerics |
| [pandas](https://pandas.pydata.org/) | 2.2.3 | BSD-3-Clause | tabular processing |
| [Matplotlib](https://matplotlib.org/) | 3.9.2 | PSF-based, BSD-compatible | every figure |
| [SciPy](https://scipy.org/) | 1.14.1 | BSD-3-Clause | statistics |
| [pysam](https://github.com/pysam-developers/pysam) | 0.22.1 | MIT | BAM iteration |
| [openpyxl](https://openpyxl.readthedocs.io/) | 3.1.5 | MIT | reading measurement spreadsheets |
| [breseq](https://github.com/barricklab/breseq) | 0.40.1 | **GPL-2.0** | contaminant divergence — see below |
| [NCBI sra-human-scrubber](https://github.com/ncbi/sra-human-scrubber) | pinned by digest | US Government work, public domain | human-read screening before deposition |
| [Kraken2](https://github.com/DerrickWood/kraken2) | via biocontainer | MIT | *optional*, prior-study reanalysis only |

## The GPL question, answered

**breseq is GPL-2.0 and this repository is MIT. That is not a conflict here.**

breseq is not vendored, not linked against, not modified, and not
redistributed. It runs as a **separate process in its own upstream
biocontainer**, pulled by digest at run time, and communicates through files on
disk. `bin/contaminant_divergence.sh` invokes it the way a person would invoke
it from a shell.

That is mere aggregation: the GPL's copyleft attaches to derivative works, and
calling a program is not deriving from it. No file in this repository is a
derivative work of breseq, so nothing here is obliged to be GPL. Anyone
redistributing the **breseq container** is bound by GPL-2.0 for that container;
this repository never does.

The same reasoning covers every other tool above, all of which are permissively
licensed anyway.

## Data redistributed by this repository

Code is not the only thing here with terms attached. Two datasets are actually
carried in the repository:

### GIAB HG002 reads — `assets/testdata/`

Redistributed under **CC0** with the NIST Data Use Policy. The full notice,
including the required source acknowledgement and the statement of how the data
was modified, is in
[`assets/testdata/README.md`](assets/testdata/README.md#redistribution-notice-for-the-giab-reads).
In short:

- source: Genome in a Bottle ultra-long ONT PromethION sequencing of
  HG002/NA24385, collected January 2019, released 8 May 2020
- **the National Institute of Standards and Technology (NIST) is acknowledged as
  the source of this data**
- **the data has been changed**: subsampled, length-filtered, partitioned by
  identity, and — for 20 records — truncated and spliced with this study's own
  reads into chimeras that do not exist in the source dataset
- consent basis: HG002/NA24385 is a Personal Genome Project participant, broadly
  consented including for commercial redistribution
- NIST makes no warranty regarding this data, and neither the modifications nor
  this repository are endorsed by NIST or UC Santa Cruz

### Prior-study values — `assets/comparison/prior_studies.tsv`

Numeric values extracted from published tables and figures, each row citing its
source publication with a DOI and the specific table or figure it came from.
These are facts reported in the literature, not copied expression, and each is
attributed at the row level so any value can be traced to its origin. The
underlying reads for Basapathi Raghavendra et al. 2023 are not redistributed
here; `bin/comparison/fetch_raghavendra.sh` retrieves them from Zenodo
(8208597).

## Reference genomes

Not redistributed. Every genome is fetched at run time by accession from NCBI
and recorded with its SHA-256 in `reference_provenance.json`, so a run is
reproducible without this repository hosting a single base of reference
sequence.
