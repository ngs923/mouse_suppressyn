# mouse_suppressyn

A pipeline for quantifying the expression of endogenous retrovirus
(ERV)-derived env genes in the mouse placenta, using public RNA-seq data and
the gEVE annotation.

---

## Contents

- [Overview](#overview)
- [Repository layout](#repository-layout)
- [Setup](#setup)
- [Running the pipeline](#running-the-pipeline)
- [Output](#output)
- [Before you run it](#before-you-run-it)
- [References](#references)

---

## Overview

```
SRA (SRR392616, SRR392617)   mouse placenta RNA-seq
  │
  ├─ Cutadapt v1.8.3         trim adapters and low-quality bases
  ├─ TopHat2 v2.1.1          align to GRCm38.p1
  ├─ Cufflinks v2.2.1        quantify with Mmus38.geve.m_v1.gtf
  └─ summarize               extract env candidates
```

## Repository layout

```
.
├── README.md
├── environment.yml           conda environment definition
├── config.example.json       template for the settings file
└── scripts/
    └── run_eve_rnaseq.py     the pipeline (Python 3)
```

Large files — the reference genome, the gEVE GTF, fastq and BAM files — are
not tracked here (see `.gitignore`). Obtain them as described below.

## Setup

### 1. Analysis tools

TopHat2 2.1.1, Cufflinks 2.2.1 and Cutadapt 1.8.3 all date from the Python 2.7
era, so keep them in a dedicated conda environment.

```bash
conda env create -f environment.yml
conda activate eve-rnaseq
```

`run_eve_rnaseq.py` itself is a Python 3 script, but the tools it calls live in
the Python 2.7 environment. Running it with that environment activated
resolves everything through `PATH`.

There are no arm64 builds of TopHat2 or Cufflinks, so on Apple Silicon Macs
build the environment for x86_64 under Rosetta, or use Docker.

```bash
CONDA_SUBDIR=osx-64 conda env create -f environment.yml
conda activate eve-rnaseq
conda config --env --set subdir osx-64
```

### 2. Reference genome

Place a GRCm38.p1 FASTA file in `ref/`. GRCm38.p1 is a 2012 patch release and
is no longer linked from the current Ensembl or NCBI landing pages; fetch it
from the Ensembl archive (around release 70) or from NCBI as
`GCA_000001635.3_GRCm38.p1`.

The choice of patch level has very little effect on the results, so the more
readily available GRCm38 primary assembly (mm10) works as a substitute.

### 3. gEVE annotation

Download `Mmus38.geve.m_v1.gtf` from the
[gEVE database](http://geve.med.u-tokai.ac.jp) and place it in `ref/`.

## Running the pipeline

```bash
cp config.example.json config.json   # edit the paths for your environment

python3 scripts/run_eve_rnaseq.py check       --config config.json
python3 scripts/run_eve_rnaseq.py inspect-gtf --config config.json
python3 scripts/run_eve_rnaseq.py all         --config config.json
```

Individual steps can be run on their own:

| Step | Description |
|---|---|
| `check` | Verify the required tools and reference files |
| `inspect-gtf` | Inspect the gEVE GTF and check chromosome names |
| `fetch` | Download fastq files from the SRA |
| `trim` | Cutadapt |
| `index` | bowtie2-build |
| `align` | TopHat2 |
| `quant` | Cufflinks |
| `summarize` | Combine FPKM values and extract env candidates |
| `all` | Run all of the above in order |

Completed steps are skipped via `.done` markers, so the same command resumes
an interrupted run. Use `--force` to redo a step and `--dry-run` to print the
commands without executing them.

## Output

```
eve_rnaseq/
├── 00_sra/             prefetched .sra files
├── 01_fastq/           fasterq-dump output (gzipped)
├── 02_trimmed/         Cutadapt output
├── 03_tophat2/<SRR>/   accepted_hits.bam and friends
├── 04_cufflinks/<SRR>/ genes.fpkm_tracking, isoforms.fpkm_tracking
├── 05_results/
│   ├── geve_fpkm_matrix.tsv     FPKM for all EVEs across all samples
│   └── geve_env_expressed.tsv   env candidates only, sorted by expression
└── logs/
```

## Before you run it

### Chromosome name mismatches

This is the most common way for the pipeline to fail. If the sequence names in
the GTF (`1` vs `chr1`) do not match the genome FASTA, **Cufflinks reports no
error and simply returns FPKM = 0 for every gene**. The `inspect-gtf` step
detects this.

```bash
python3 scripts/run_eve_rnaseq.py inspect-gtf --config config.json --fix-chrom-names
```

With `--fix-chrom-names`, a converted copy of the GTF matching the FASTA is
written out.

### Multi-mapping reads

ERVs are repetitive, so this choice strongly affects the results. The defaults
are TopHat2 `-g 20` (up to 20 locations per read) and Cufflinks `-u`
(multi-read rescue).

If you intend to make claims about a specific env locus, also check the result
using uniquely mapping reads only (`-g 1`). Both behaviours are controlled by
`tophat_max_multihits` and `cufflinks_multi_read_correct`.

### Adapter sequences

The defaults are Illumina TruSeq adapters. SRR392616 and SRR392617 are from
2012, so the library preparation of the day may have used different adapters.
Run FastQC on the raw fastq files and check the adapter content before
settling on these values.

### Single-end vs paired-end

The layout is detected from the number of files produced by
`fasterq-dump --split-3`, and the Cutadapt and TopHat2 arguments are adjusted
accordingly. For paired-end data, tune `tophat_mate_inner_dist` to match your
library (the default of 50 is a placeholder).

### The env pattern

`summarize` treats any gene whose `gene_id` contains `env` as an env
candidate. Check the gEVE naming convention in the `inspect-gtf` output and
adjust `env_pattern` if needed.

## References

- Nakagawa S, Takahashi MU. gEVE: a genome-based endogenous viral element database provides comprehensive viral protein-coding sequences in mammalian genomes. *Database* 2016; baw087.
- Kim D, Pertea G, Trapnell C, et al. TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions. *Genome Biol* 2013; 14: R36.
- Trapnell C, Williams BA, Pertea G, et al. Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation. *Nat Biotechnol* 2010; 28: 511–515.
- Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet J* 2011; 17: 10–12.
