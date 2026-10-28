# Spatial transcriptomics

In this practical, we explore how spatial transcriptomics can reveal the organization of gene expression within mouse liver tissue.

Unlike conventional single-cell or single-nucleus RNA-seq, spatial transcriptomics preserves information about the anatomical location of transcripts.

We compare three spatial transcriptomics technologies using publicly available mouse liver datasets.

## Datasets

### 10x Visium

> Guilliams M, et al. *Spatial proteogenomics reveals distinct and evolutionarily conserved hepatic macrophage niches.* Cell (2022).
> https://doi.org/10.1016/j.cell.2021.12.018

We investigate spatial gene-expression patterns using sequencing-based transcriptomics and histological images.

Data: [GSE192741](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192741)

### MERFISH

> Liu J, et al. *Concordance of MERFISH spatial transcriptomics with bulk and single-cell RNA sequencing.* Life Science Alliance (2022).
> https://doi.org/10.26508/lsa.202201701

We explore imaging-based spatial transcriptomics, cell segmentation, clustering and spatial gene-expression analysis.

### Resolve Molecular Cartography

> Bravo González-Blas C, et al. *Single-cell spatial multi-omics and deep learning dissect enhancer-driven gene regulatory networks in liver zonation.* Nature Cell Biology (2024).
> https://doi.org/10.1038/s41556-023-01316-4

We investigate cell identities and liver zonation using a targeted imaging-based spatial transcriptomics dataset.

## Notebook

[Spatial transcriptomics analysis of the mouse liver](Regulatory_genomics_example_spatial_transcriptomics_of_the_mouse_liver.ipynb)

The practical covers:

1. Loading and exploring spatial transcriptomics datasets.
2. Quality control and normalization.
3. Dimensionality reduction and clustering.
4. Visualizing gene expression in tissue.
5. Image segmentation using thresholding and watershed.
6. Identifying spatially variable genes using Moran's I.
7. Investigating spatial neighborhoods and liver zonation.

## Software

The practical uses Python, primarily Scanpy, Squidpy, AnnData, NumPy, SciPy and scikit-image.

Instructions for creating the `regulatory_genomics` Conda environment are included in the notebook.

The same environment is used for all three spatial transcriptomics technologies.

## Data

The Visium and MERFISH datasets are downloaded automatically by the notebook.

The Resolve Molecular Cartography dataset must be downloaded manually from [SCope](https://scope.aertslab.org/#/Bravo_et_al_Liver).

Instructions are provided in the notebook.
