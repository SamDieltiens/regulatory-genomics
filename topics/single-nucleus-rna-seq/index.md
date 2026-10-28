# Single-nucleus RNA-seq

In this practical, we analyze a published single-nucleus RNA-seq dataset from the mouse liver:

> Bravo González-Blas C, et al. *Single-cell spatial multi-omics and deep learning dissect enhancer-driven gene regulatory networks in liver zonation.* Nature Cell Biology (2024).
> https://doi.org/10.1038/s41556-023-01316-4

The liver contains different cell populations, including hepatocytes, endothelial cells, Kupffer cells and hepatic stellate cells. Hepatocytes also exhibit transcriptional heterogeneity associated with their location along the portal–central axis, a phenomenon known as liver zonation.

In this practical, we investigate whether cell identities and hepatocyte zonation-associated transcriptional programs can be recovered from single-nucleus RNA-seq data.

## Notebook

[Single-nucleus RNA-seq analysis of the mouse liver](Regulatory_genomics_example_single_nucleus_RNA_seq_of_the_mouse_liver.ipynb)

The analysis covers:

1. Downloading and exploring the gene-expression count matrix.
2. Quality control and filtering.
3. Normalization and highly variable gene selection.
4. Dimensionality reduction using PCA and UMAP.
5. Leiden clustering and cell-type annotation.
6. Investigating hepatocyte heterogeneity and liver zonation.

## Software

The analysis uses Python, primarily Scanpy, AnnData, NumPy, Pandas and Matplotlib.

The notebook contains instructions for creating the `regulatory_genomics` Conda environment and registering the corresponding Jupyter kernel.

## Data

The published single-nucleus RNA-seq data are available from GEO:

[GSE218465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218465)

The notebook automatically downloads the processed count matrix. No manual data download is required.
