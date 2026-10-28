# MPRA (massively parallel reporter assay)

In this practical we reproduce the MPRA analysis of:

> Bravo González-Blas C, Matetovici I, Hillen H, et al. *Single-cell spatial multi-omics and deep learning dissect enhancer-driven gene regulatory networks in liver zonation.* Nature Cell Biology (2024) 26:153–167. https://www.nature.com/articles/s41556-023-01316-4

The authors tested thousands of candidate enhancers of the mouse liver in a reporter assay, both in the living mouse liver and in the human liver cell line HepG2, to find out which accessible regions are really active enhancers and which transcription factor motifs make them active. The data is available on GEO: [GSE218472](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218472).

The practical consists of two notebooks:

1. [Part 1: from sequencing reads to barcode counts](MPRA_part_1_reads_to_counts.ipynb) (bash): downloading the data (SRA and ENA), quality control (FastQC, MultiQC), cutting barcodes out of the reads (cutadapt, fastp, seqkit), building the enhancer-barcode assignment and counting barcodes per enhancer.
2. [Part 2: enhancer activity and biological interpretation](MPRA_part_2_activity_and_interpretation.ipynb) (Python): quality control of the count tables, modelling cDNA versus plasmid counts with DESeq2 (PyDESeq2), calling active enhancers against shuffled negative controls, comparing activity in vivo and in HepG2, motif enrichment (HOMER) and the effect of mutating transcription factor binding sites, including enhancer zonation measured on FACS-sorted hepatocytes.

## Software

All software for both notebooks is installed in a single mamba environment. The exact commands (including how to register the environment as a Jupyter kernel) are given at the start of each notebook.

    mamba create -n mpra_teaching -c conda-forge -c bioconda python=3.11 bash_kernel ipykernel \
        sra-tools fastqc multiqc cutadapt fastp seqkit pigz \
        pandas numpy scipy statsmodels matplotlib seaborn pydeseq2 openpyxl homer bedtools

| Notebook | Kernel |
|---|---|
| Part 1 | Bash (liver MPRA) |
| Part 2 | Python (liver MPRA) |

The motif analysis in part 2 needs the mouse genome sequence (mm10). On the VSC cluster a copy is available in the shared genome folder (`/lustre1/project/res_00001/genomes/`); otherwise download it once from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz).

## Data

All data is stored in the folder `data`, next to the notebooks. Part 1 downloads ~ 5 GB of raw sequencing reads; part 2 only needs the (small) count tables from GEO and the supplementary tables of the paper, so **part 2 can be run without running part 1**.
