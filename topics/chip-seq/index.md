# ChIP-seq

In this practical we analyse the CEBPA ChIP-seq data of:

> Ballester B, Medina-Rivera A, Schmidt D, et al. *Multi-species, multi-transcription factor binding highlights conserved control of tissue-specific biological pathways.* eLife (2014) 3:e02626. https://elifesciences.org/articles/02626

The authors mapped the binding sites of four liver transcription factors (CEBPA, HNF4A, FOXA1 and ONECUT1/HNF6) in the livers of five species using ChIP-seq. We use the human CEBPA ChIP-seq experiment and its control (Input) to find out where CEBPA binds, which DNA motif it recognises and which genes and processes it controls. The data is available on ENA: [PRJEB1571](https://www.ebi.ac.uk/ena/browser/view/PRJEB1571) (ArrayExpress E-MTAB-1509).

To find out whether the CEBPA bound regions are active enhancers, we combine them with H3K27ac ChIP-seq data of liver and spleen of the same donor from the [ENTEx](https://www.encodeproject.org/search/?type=Experiment&internal_tags=ENTEx) project (ENCODE: ENCSR864OOO, ENCSR343DFX and their Inputs ENCSR309LQG, ENCSR749YQP).

The practical consists of one notebook:

1. [ChIP-seq: where does the transcription factor CEBPA bind in the human liver?](ChIP-seq.ipynb) (bash): downloading the data (ENA), quality control (FastQC, MultiQC), aligning reads to the genome (Bowtie2), filtering alignments and removing duplicates (samtools), visualisation in IGV and coverage tracks (deepTools), ChIP quality (fingerprint plot, FRiP), peak calling (MACS3), peak annotation, motif analysis and gene ontology enrichment (HOMER), and the comparison with H3K27ac in liver and spleen, including a differential analysis (MACS3 bdgdiff).

## Software

All software is installed in a single mamba environment. The exact commands (including how to register the environment as a Jupyter kernel) are given at the start of the notebook.

    mamba create -n chipseq -c conda-forge -c bioconda --strict-channel-priority \
        fastqc=0.12.1 multiqc=1.35 bowtie2=2.5.4 samtools=1.22.1 bedtools=2.31.1 \
        deeptools=3.5.6 macs3=3.1.0 homer=5.1 bash_kernel curl bc

HOMER also needs its copy of the human genome (hg38, ~ 3 GB), which is installed once with `perl ${CONDA_PREFIX}/share/homer/configureHomer.pl -install hg38`.

| Notebook | Environment | Kernel |
|---|---|---|
| ChIP-seq | `chipseq` | Bash (ChIP-seq) |

The reads are aligned to hg38 with a Bowtie2 index. On the VSC cluster an index is available at `/lustre1/project/res_00001/genomes/homo_sapiens/hg38_iGenomes/indexes/Bowtie2/genome`; elsewhere it can be built from the genome sequence with `bowtie2-build` (see Step 3 of the notebook).

## Data

All data is stored in the folder `data`, next to the notebook. Running the full analysis requires ~ 30 GB of disk space (most of it for the H3K27ac alignments and the bedGraph files of MACS3). The full notebook runs in ~ 45 minutes using 8 threads, including ~ 15 minutes to download the data.
