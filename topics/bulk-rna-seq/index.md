# Bulk RNA-seq

In this practical we reproduce the bulk RNA-seq analysis of:

> Hansen JL, Loell KJ, Cohen BA. *A test of the pioneer factor hypothesis using ectopic liver gene activation.* eLife (2022) 11:e73358. https://elifesciences.org/articles/73358

The authors expressed the transcription factors FOXA1, HNF4A, or both, in K562 cells (a blood cell line) and measured gene expression using RNA-seq. The data is available on GEO: [GSE182190](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182190).

The practical consists of two notebooks:

1. [Part 1: from sequencing reads to a count matrix](part_1_reads_to_counts.ipynb) (bash): downloading the data, quality control (FastQC, MultiQC), aligning reads to the genome (STAR), working with SAM/BAM files (samtools) and counting reads per gene (featureCounts).
2. [Part 2: differential expression and biological interpretation](part_2_differential_expression.ipynb) (R): exploratory analysis (PCA), differential expression (DESeq2), tissue enrichment using the Human Protein Atlas and gene ontology enrichment (clusterProfiler).

## Software

The software for each notebook is installed in a separate mamba environment, the exact commands (including how to register the environment as a Jupyter kernel) are given at the start of each notebook.

| Notebook | Environment | Kernel |
|---|---|---|
| Part 1 | `rnaseq_tools` | Bash (FoxA1/Hnf4a RNA-seq) |
| Part 2 | `rnaseq_r` | R (FoxA1/Hnf4a RNA-seq) |

## Data

All data is stored in the folder `data`, next to the notebooks. Running the full analysis requires ~ 200 GB of disk space and ~ 32 GB of memory (to build the STAR index).
