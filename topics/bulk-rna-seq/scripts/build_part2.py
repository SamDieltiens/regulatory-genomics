"""Build (and optionally execute) part_2_differential_expression.ipynb"""
import sys
import os
import nbformat as nbf
TOPIC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(TOPIC)), "tools"))
from nbtools import md, code, new_notebook, execute, sync_markdown

OUT = os.path.join(TOPIC, "part_2_differential_expression.ipynb")

cells = [
md("""
# Bulk RNA-seq part 2: differential expression and biological interpretation

**Case study: Can the transcription factors FOXA1 and HNF4A activate liver genes in a blood cell line?**

Hansen JL, Loell KJ, Cohen BA. *A test of the pioneer factor hypothesis using ectopic liver gene activation.* eLife (2022) 11:e73358. https://elifesciences.org/articles/73358

K562 cells (a blood cell line) with doxycycline inducible expression of FOXA1, HNF4A or both FOXA1 and HNF4A. For each line we have 3 uninduced and 3 induced (0.5 µg/ml doxycycline, 24 hours) samples.

**Research questions:**

1. Which genes are activated upon expression of FOXA1, HNF4A or both?
2. Are these genes liver (or intestine) specific?
3. Do FOXA1 and HNF4A activate the same genes or different genes?

Here we will use DESeq2 to perform differential gene expression analysis. DESeq2 has a nice tutorial that can be found [here](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

**Note**:

All R packages can be installed in a single mamba environment. This environment also contains `r-irkernel`, which allows us to run R code in a Jupyter notebook (and `jupyter_client`, which is needed to register the kernel).

    mamba create -n rnaseq_r -c conda-forge -c bioconda r-base=4.4 bioconductor-deseq2 bioconductor-apeglm bioconductor-clusterprofiler bioconductor-org.hs.eg.db bioconductor-enrichplot r-ggplot2 r-ggrepel r-pheatmap r-dplyr r-tidyr r-irkernel r-upsetr jupyter_client

To use this environment in Jupyter we have to register it as a *kernel*. Run the following commands in a terminal (not in this notebook):

```
mamba activate rnaseq_r
R -e 'IRkernel::installspec(name = "ir-foxa1_hnf4a", displayname = "R (FoxA1/Hnf4a RNA-seq)")'
jupyter kernelspec list   # the kernel ir-foxa1_hnf4a should be listed
```

Next, (re)start Jupyter and select the kernel "R (FoxA1/Hnf4a RNA-seq)" (Kernel > Change Kernel).
"""),
md("""
## Step 0. set up working directory

We will use the count matrix generated in the previous notebook, which is stored in the `data` folder.
"""),
code("""
setwd("data")
options(repr.plot.width = 7, repr.plot.height = 5) # size of the plots in the notebook
"""),
code("""
suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(pheatmap)
})
"""),
md("""
## Step 1. Load data

The starting point of a DESeq2 analysis is a count matrix K with one row for each gene i and one column for each sample j. The matrix entries $K_{ij}$ indicate the number of sequencing reads (here read pairs) that have been unambiguously mapped to a gene in a sample. Let's load the count matrix we generated in the previous notebook.
"""),
code("""
cts <- as.matrix(read.csv("counts/foxa1_hnf4a.counts.formatted", sep="\\t", row.names = "Geneid"))
"""),
code("""
head(cts, 5)
"""),
code("""
dim(cts)
"""),
md("""
We have 18 samples and more than 60,000 genes (this includes protein coding genes, but also non-coding RNAs, pseudogenes, ...).

Next we will generate a dataframe containing metadata on our samples. We can extract this information from the sample names (e.g. `FOXA1HNF4A_I2` is the FOXA1-HNF4A line, induced, replicate 2).
"""),
code("""
coldata <- data.frame(
    row.names = colnames(cts),
    line      = sub("_.*", "", colnames(cts)),                                       # everything before the "_"
    condition = ifelse(grepl("_U", colnames(cts)), "uninduced", "induced"),          # U = uninduced, I = induced
    replicate = sub(".*_[UI]", "rep", colnames(cts))                                 # the number after U or I
)
coldata$line      <- factor(coldata$line, levels = c("FOXA1", "HNF4A", "FOXA1HNF4A"))
coldata$condition <- factor(coldata$condition, levels = c("uninduced", "induced"))  # uninduced is our baseline
coldata$group     <- factor(paste(coldata$line, coldata$condition, sep = "_"))
"""),
code("""
coldata
"""),
md("""
**It's critical that the columns of the count matrix and the rows of the coldata are in the same order**. This is the case here.
"""),
code("""
all(rownames(coldata) == colnames(cts))
"""),
md("""
## Step 2. Exploratory analysis: how similar are the samples?

Before testing for differential expression it is a good idea to check how similar the samples are to each other. We expect that replicates are very similar, and that induced and uninduced samples are different.

First we construct a DESeqDataSet containing all 18 samples. For this we need:
- the count matrix (cts)
- the metadata (coldata)
- the design, this specifies which variables to use to perform differential gene expression testing. Here we use the `group` column (the combination of cell line and condition).
"""),
code("""
dds_all <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~ group)
dds_all
"""),
md("""
We will perform some prefiltering to only keep genes that have a count of at least 10 for a minimum of 3 samples (the size of our smallest group).
"""),
code("""
smallestGroupSize <- 3
keep <- rowSums(counts(dds_all) >= 10) >= smallestGroupSize
dds_all <- dds_all[keep,]
dds_all
"""),
md("""
Raw counts are not suitable to compare samples:

- samples are sequenced to a different depth (the FOXA1-HNF4A samples have twice as many reads!)
- the variance of counts depends on the mean: highly expressed genes have a much larger variance than lowly expressed genes. If we would compute distances on raw counts, they would be dominated by a few highly expressed genes.

DESeq2 provides the *variance stabilizing transformation* (`vst`). This normalizes the counts for sequencing depth and transforms them to a log2-like scale where the variance is (approximately) independent of the mean.
"""),
code("""
vsd_all <- vst(dds_all, blind = TRUE)
head(assay(vsd_all), 3)
"""),
md("""
### Principal component analysis (PCA)

A PCA summarizes the expression of thousands of genes into a few *principal components* (PCs). The first PC captures the most variation between the samples, the second PC captures the most variation that is not captured by the first PC, etc. Samples that are close to each other in a PCA plot have a similar expression profile.

DESeq2 has a function `plotPCA` which by default uses the 500 most variable genes.
"""),
code("""
pcaData <- plotPCA(vsd_all, intgroup = c("line", "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = line, shape = condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_bw()
"""),
md("""
Each dot is a sample. Note that:

- Replicates are very close to each other (the three dots of each group almost overlap), this means that the replicates are very reproducible.
- The uninduced samples of the three lines are close to each other (right side of the plot): without doxycycline the three lines are very similar, as expected (they are all K562 cells).
- The first principal component (PC1, 63% of the variance) separates the induced samples of the HNF4A and FOXA1-HNF4A lines from the uninduced samples. The induced FOXA1 samples are much closer to the uninduced samples on PC1.
- The second principal component (PC2) separates the induced samples of the three lines from each other: FOXA1, HNF4A and FOXA1-HNF4A each induce a different expression program.

**Question**: Based on this plot, which transcription factor has the largest effect on gene expression, FOXA1 or HNF4A?
"""),
md("""
### Sample-to-sample distances

Another way to look at this is by calculating the (euclidean) distance between all pairs of samples and plotting these as a heatmap.
"""),
code("""
sampleDists <- dist(t(assay(vsd_all)))
pheatmap(as.matrix(sampleDists),
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         color = colorRampPalette(c("navy", "white"))(100),
         annotation_col = coldata[, c("line", "condition")])
"""),
md("""
The samples first cluster by condition (induced vs uninduced) and then by cell line. Again we see that replicates are very similar (dark blue squares along the diagonal).

Because the three cell lines are separate clonal lines with a different baseline expression, we will compare induced vs uninduced samples **within each cell line**. We will start with the FOXA1-HNF4A line (Figure 1 of the paper) and later repeat the analysis for the FOXA1 and HNF4A lines (Figure 2 of the paper).
"""),
md("""
## Step 3. Model read counts and performing differential expression analysis

Let's make a DESeqDataSet containing only the FOXA1-HNF4A samples. Our design is now `~ condition`, as we want to test for differential expression between induced and uninduced.
"""),
code("""
samples_fh <- rownames(coldata)[coldata$line == "FOXA1HNF4A"]
dds <- DESeqDataSetFromMatrix(countData = cts[, samples_fh],
                              colData = coldata[samples_fh, ],
                              design = ~ condition)
dds$condition <- factor(dds$condition, levels = c("uninduced","induced")) # uninduced is our baseline
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds
"""),
md("""
Given the count matrix: $K_{ij}$ with one row for each gene $i$ and one column for each sample $j$

Read counts $K_{ij}$ are modeled as a negative binomial distribution (a.k.a. gamma-Poisson distribution) with a fitted mean $\\mu_{ij}$ and gene specific dispersion $\\alpha_{i}$.

<p style="text-align: center;">$K_{ij} \\sim \\text{NB}(\\mu_{ij}, \\alpha_{i})$</p>

For illustrative purposes, the distribution of the expression of genes in sample FOXA1HNF4A_U1 is shown.
"""),
code("""
hist(counts(dds)[, 'FOXA1HNF4A_U1'], breaks = 1000, xlim = c(0, 5000), col = "darkgray", main = "FOXA1HNF4A_U1", xlab = "counts")
"""),
md("""
Most genes have a low count, while a few genes have a very high count.

The fitted mean $\\mu_{ij}$ is composed of sample-specific size factor $s_j$ and a parameter $q_{ij}$, proportional to the concentration of cDNA fragments from the gene $i$ in sample $j$.

<p style="text-align: center;">$\\mu_{ij}=s_j \\cdot q_{ij}$</p>

The parameter $s_j$, which is the same for all genes in a sample, accounts for differences in sequencing depths between two samples. It can also be advantageous to calculate a gene-specific normalization factor $s_{ij}$. This gene- and sample-specific normalization factor accounts for further sources of technical bias such as differing dependence on GC content, gene length, ...

In this case the formula looks like this:

<p style="text-align: center;">$\\mu_{ij}=s_{ij} \\cdot q_{ij}$</p>

Finally generalized linear models (GLMs) are used to fit for each gene as follows:

<p style="text-align: center;">$\\log_2{q_{ij}} = \\sum_r{x_{jr} \\cdot \\beta_{ir} }$</p>

With $x_{jr}$ a design matrix for $r$ conditions indicating whether sample $j$ belongs to each condition $r$ or not.

In this case we have only a single condition (uninduced vs induced) and thus the formula simplifies to:

<p style="text-align: center;">$\\log_2{q_{ij}} = \\beta_{i0} + x_{j} \\cdot \\beta_{i} $</p>

with $x_j = 1$ for induced samples and $x_j = 0$ for uninduced samples. $\\beta_{i0}$ is the (log2) expression in uninduced samples and $\\beta_{i}$ gives the log2 fold change for gene $i$.
"""),
md("""
The variability between replicates is modeled using the dispersion parameter $\\alpha_{i}$. It defines the relationship between the variance of the observed count and its mean value (how far do we expect the observed count will be from the mean value).

<p style="text-align: center;">$\\text{Var}(K_{ij}) = E[(K_{ij} - \\mu_{ij})^2] = \\mu_{ij} + \\alpha_i \\cdot \\mu_{ij}^2$</p>

Accurate estimation of this dispersion parameter $\\alpha_i$ is critical for statistical testing of differential expression. For studies with a small sample size this parameter is difficult to estimate. DESeq2 solves this issue by assuming that genes with similar average expression values will have a similar dispersion.

Let's plot the mean vs the variance of the counts of the uninduced samples.
"""),
code("""
cts_u <- counts(dds)[, dds$condition == "uninduced"]
data <- data.frame(
    mean = rowMeans(cts_u),
    var  = rowVars(cts_u))
ggplot(data, aes(x = mean, y = var)) +
    geom_point(size = 0.5) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 2)
"""),
md("""
The Red line represents the case where the variance equals the mean (this would be the case if the counts followed a Poisson distribution). For highly expressed genes the variance is much larger than the mean, this is called *overdispersion*.
"""),
md("""
The **size factors** $s_{j}$ can be estimated for each sample using the function

        estimateSizeFactors
"""),
code("""
dds <- estimateSizeFactors(dds)
"""),
code("""
dds@colData
"""),
md("""
The size factors are roughly proportional to the sequencing depth of each sample.
"""),
code("""
data.frame(sizeFactor = sizeFactors(dds), total_counts = colSums(counts(dds)))
"""),
md("""
The **dispersion** $\\alpha_{i}$ can be estimated using the function

    estimateDispersions
"""),
code("""
dds <- estimateDispersions(dds)
"""),
code("""
head(mcols(dds))
"""),
code("""
ggplot(as.data.frame(mcols(dds)[!mcols(dds)$allZero, ])) +
    geom_point(aes(x = baseMean, y = baseVar), size = 0.5) +
    geom_line(aes(x = baseMean, y = baseMean + dispFit * baseMean * baseMean ), color = "red", linewidth = 2) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10")
"""),
md("""
As you can see the red line now fits the data much better!

DESeq2 also has a built in function to plot the dispersion estimates: the black dots are the gene-wise estimates, the red line is the fitted trend and the blue dots are the final estimates, which are *shrunk* towards the trend line.
"""),
code("""
plotDispEsts(dds)
"""),
md("""
Finally we can **fit the GLM for each gene** to obtain $\\beta_i$, the log2FC of each gene across the two conditions, using the function:

    nbinomWaldTest
"""),
code("""
dds <- nbinomWaldTest(dds)
"""),
md("""
The results of this can be visualized using the function

    results
"""),
code("""
res <- results(dds)
res
"""),
md("""
- baseMean: the average of the normalized counts over all samples
- log2FoldChange: the log2 fold change (induced vs uninduced, see the header of the output)
- lfcSE: the standard error of the log2 fold change
- stat: the Wald statistic (log2FoldChange / lfcSE)
- pvalue: the p-value of the Wald test
- padj: the p-value corrected for multiple testing (Benjamini-Hochberg). We test more than 15,000 genes, so by chance alone 5% of them (> 750 genes!) would have a p-value < 0.05.

The three functions above:

    estimateSizeFactors
    estimateDispersions
    nbinomWaldTest

can be run in one go using the wrapper function:

    DESeq
"""),
code("""
dds <- DESeq(dds)
res <- results(dds)
res
"""),
code("""
summary(res, alpha = 0.05)
"""),
md("""
Genes with a low read count have a strong variance of the log2FoldChange (see plot below). Because of this, weakly expressed genes can seem to have strong differences between two conditions. However, this is not due to a real difference but due to noise in the measurements. This will cause issues if we want to rank genes based on their log2fold change values.
"""),
code("""
ggplot(as.data.frame(res), aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(colour = padj < 0.05), size = 0.5) +
    scale_x_continuous(trans = "log10") +
    scale_color_manual(name = "padj < 0.05", values = setNames(c("red", "black"), c(T, F)))
"""),
md("""
As you can see genes with a low mean (left of the plot) have a big variance of the log2foldchange. This can also be seen by looking at the standard error of the log2foldchange vs the mean
"""),
code("""
ggplot(as.data.frame(res), aes(x = baseMean, y = lfcSE)) +
    geom_point(size = 0.5) +
    scale_x_continuous(trans = "log10") +
    ylim(0, 3)
"""),
md("""
To account for this DESeq2 has a function to *shrink* Log2FoldChange towards zero in such a manner that the shrinkage is stronger when the available information for a gene is low (i.e. low number of counts, or dispersion is high).

The function for this is called

    lfcShrink

it takes the `dds` object as input together with the coefficient we want to shrink
"""),
code("""
resultsNames(dds)
"""),
md("""
We want to shrink the coefficient `condition_induced_vs_uninduced` as this stores our Log2FoldChange values.
"""),
code("""
resLFC <- lfcShrink(dds, coef="condition_induced_vs_uninduced", type = "apeglm")
"""),
code("""
resLFC
"""),
code("""
ggplot(as.data.frame(resLFC), aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(colour = padj < 0.05), size = 0.5) +
    scale_x_continuous(trans = "log10") +
    scale_color_manual(name = "padj < 0.05", values = setNames(c("red", "black"), c(T, F)))
"""),
md("""
As you can see, now the effect of the mean on the variance of the Log2FoldChange is much smaller. Again we can also visualize this by plotting the standard error.
"""),
code("""
ggplot(as.data.frame(resLFC), aes(x = baseMean, y = lfcSE)) +
    geom_point(size = 0.5) +
    scale_x_continuous(trans = "log10") +
    ylim(0, 3)
"""),
md("""
## Step 4. Downstream analysis
"""),
md("""
### Are the transgenes induced?

Before we look at other genes, let's check whether FOXA1 and HNF4A are indeed induced upon doxycycline treatment. The function `plotCounts` plots the normalized counts of a single gene.
"""),
code("""
resLFC_df <- as.data.frame(resLFC)
resLFC_df[c("FOXA1", "HNF4A"), ]
"""),
code("""
options(repr.plot.width = 10, repr.plot.height = 4)
d <- rbind(
    cbind(plotCounts(dds, gene = "FOXA1", intgroup = "condition", returnData = TRUE), gene = "FOXA1"),
    cbind(plotCounts(dds, gene = "HNF4A", intgroup = "condition", returnData = TRUE), gene = "HNF4A"))
ggplot(d, aes(x = condition, y = count)) +
    geom_point(size = 3) +
    scale_y_log10() +
    facet_wrap(~ gene) +
    ylab("normalized counts") +
    theme_bw()
options(repr.plot.width = 7, repr.plot.height = 5)
"""),
md("""
HNF4A is very strongly induced (log2 fold change of 9.5, i.e. ~ 700 fold). As we have seen in the previous notebook, we can't measure the expression of the FOXA1 transgene because it is codon optimized: the reads coming from the transgene do not map to the *FOXA1* gene. The reads that map to *FOXA1* come from the endogenous gene, which is even downregulated upon induction.
"""),
md("""
### Getting a list of up- or downregulated genes
To manipulate dataframes we use the package dplyr.

Let's define some thresholds
"""),
code("""
log2FC_thr = 1
pad_thr = 0.05
"""),
code("""
upregulated_genes <- row.names(resLFC_df %>%
    filter(log2FoldChange > log2FC_thr, padj < pad_thr))
downregulated_genes <- row.names(resLFC_df %>%
    filter(log2FoldChange < -log2FC_thr, padj < pad_thr))
length(upregulated_genes)
length(downregulated_genes)
"""),
md("""
That's a lot of genes! In fact, without a log2 fold change threshold more than 75% of the genes have a padj < 0.05 (see the output of `summary(res)` above). This is because these samples are sequenced very deeply and the replicates are very similar (low dispersion), which gives us a lot of statistical power: even very small changes in expression become significant. A significant change is not necessarily a biologically relevant change. This is why we also use a threshold on the log2 fold change.

### Sorting the gene list based on log2FC values
"""),
code("""
resLFC_df %>%
    filter(padj < pad_thr) %>%
    arrange(desc(log2FoldChange)) %>%
    head(20)
"""),
md("""
Many of the genes with the highest log2 fold change are typical liver genes, for example the apolipoproteins APOA1, APOB and APOC3, which are made by the liver (and the intestine) and transport lipids through the blood. Also HNF1A, another liver transcription factor, is activated.

Note that for many of these genes the lfcSE (the uncertainty on the log2 fold change) is large (~3). These genes have (almost) zero counts in the uninduced samples, so the fold change is very uncertain: is it 2^12 or 2^15 fold? It doesn't really matter, the gene is clearly switched on!

### Volcano plot

A volcano plot shows the log2 fold change (x-axis) versus the significance (-log10 adjusted p-value, y-axis) of each gene.
"""),
code("""
resLFC_df$gene <- rownames(resLFC_df)
resLFC_df$status <- "not significant"
resLFC_df$status[resLFC_df$gene %in% upregulated_genes]   <- "up"
resLFC_df$status[resLFC_df$gene %in% downregulated_genes] <- "down"
# some padj values are so small that they are rounded to 0 (-log10(0) is infinite), we set these to 1e-300
resLFC_df$padj_plot <- pmax(resLFC_df$padj, 1e-300)
genes_to_label <- c("HNF4A", "FOXA1", "ALB", "APOA1", "APOB", "APOC3", "SERPINA1", "HNF1A", "TTR", "AFP")
top_genes <- resLFC_df %>% filter(gene %in% genes_to_label)
ggplot(resLFC_df %>% filter(!is.na(padj)), aes(x = log2FoldChange, y = -log10(padj_plot), color = status)) +
    geom_point(size = 0.5) +
    geom_text_repel(data = top_genes, aes(label = gene), color = "black", size = 3, max.overlaps = 20) +
    scale_color_manual(values = c("up" = "red", "down" = "blue", "not significant" = "gray")) +
    geom_vline(xintercept = c(-log2FC_thr, log2FC_thr), linetype = "dashed") +
    geom_hline(yintercept = -log10(pad_thr), linetype = "dashed") +
    ylab("-log10(padj)") +
    theme_bw()
"""),
md("""
Many genes are strongly upregulated, including several well known liver genes (ALB, APOA1, APOB, SERPINA1, ...). Note that more genes are strongly *up* than strongly *down* regulated (the plot is asymmetric): FOXA1 and HNF4A are transcriptional activators.

### Plotting a heatmap

To plot a heatmap we use the variance stabilized counts (see above) and scale each gene (row) to have a mean of 0 and a standard deviation of 1 (`scale = "row"`). This way we can compare the pattern of genes with very different expression levels.
"""),
code("""
vsd <- vst(dds, blind = FALSE)
pheatmap(assay(vsd)[c(upregulated_genes, downregulated_genes), ],
         scale = "row", show_rownames = F,
         annotation_col = coldata[, c("condition"), drop = FALSE])
"""),
md("""
The heatmap shows the expression of all up- and downregulated genes: the replicates are very consistent.
"""),
md("""
## Step 5. Biological interpretation

### Which genes are *activated*?

We are mostly interested in genes that are silent in K562 cells and that become activated by FOXA1 and HNF4A. In the paper (see Methods) a gene was called:

- *upregulated* if it had a log2 fold change of at least 1
- *activated* if, in addition, it had fewer than 50 normalized reads in the uninduced control.

We will use the same definition, but we will also require that the change is significant (padj < 0.05).

Let's write a function that gives us the activated genes, so we can easily reuse this later for the other cell lines. Normalized counts can be obtained with `counts(dds, normalized = TRUE)`.
"""),
code("""
get_activated_genes <- function(dds, resLFC, log2FC_thr = 1, pad_thr = 0.05, max_uninduced = 50) {
    norm_counts <- counts(dds, normalized = TRUE)
    mean_uninduced <- rowMeans(norm_counts[, dds$condition == "uninduced"])
    df <- as.data.frame(resLFC)
    df$mean_uninduced <- mean_uninduced[rownames(df)]
    rownames(df %>% filter(log2FoldChange >= log2FC_thr, padj < pad_thr, mean_uninduced < max_uninduced))
}
activated_fh <- get_activated_genes(dds, resLFC)
length(activated_fh)
"""),
md("""
More than 2000 genes are activated: these genes are (almost) silent in uninduced K562 cells, and are switched on by FOXA1 and HNF4A.
"""),
md("""
### Are the activated genes liver specific?

FOXA1 and HNF4A are important TFs in the liver (and intestine). Are the activated genes liver (or intestine) specific genes?

To answer this question we need a list of tissue specific genes. The [Human Protein Atlas](https://www.proteinatlas.org/) (HPA) measured the expression of all genes in many human tissues using RNA-seq. A gene is called *tissue enriched* if it has at least four-fold higher expression in one tissue compared to all other tissues. This is the same definition as used in the paper.

Let's download the HPA data.
"""),
code("""
if (!file.exists("annotation/proteinatlas.tsv")) {
    dir.create("annotation", showWarnings = FALSE)
    download.file("https://www.proteinatlas.org/download/proteinatlas.tsv.zip", "annotation/proteinatlas.tsv.zip")
    unzip("annotation/proteinatlas.tsv.zip", exdir = "annotation")
}
hpa <- read.csv("annotation/proteinatlas.tsv", sep = "\\t", check.names = FALSE)
dim(hpa)
"""),
code("""
hpa[hpa$Gene %in% c("ALB", "FOXA1", "ACTB"), c("Gene", "Ensembl", "RNA tissue specificity", "RNA tissue specific nTPM")]
"""),
md("""
- `RNA tissue specificity`: the category of tissue specificity (Tissue enriched, Group enriched, Tissue enhanced, Low tissue specificity, Not detected)
- `RNA tissue specific nTPM`: the tissue(s) in which the gene is specifically expressed and the expression level (nTPM: normalized transcripts per million)

ALB is liver enriched, FOXA1 is higher expressed in the prostate than in other tissues (tissue enhanced) and ACTB (beta-actin) is expressed in all tissues (low tissue specificity).

Let's get the tissue enriched genes for the liver and the intestine.
"""),
code("""
hpa_enriched <- hpa[hpa$`RNA tissue specificity` == "Tissue enriched", ]
hpa_enriched$tissue <- sub(":.*", "", hpa_enriched$`RNA tissue specific nTPM`)  # remove the expression value
head(sort(table(hpa_enriched$tissue), decreasing = TRUE), 10)
"""),
code("""
liver_genes     <- hpa_enriched$Gene[hpa_enriched$tissue == "liver"]
intestine_genes <- hpa_enriched$Gene[hpa_enriched$tissue == "intestine"]
length(liver_genes)
length(intestine_genes)
"""),
md("""
The current version of the HPA has 267 liver enriched genes and 123 intestine enriched genes (the paper used an older version with 242 and 122 genes).

How many of the activated genes are liver or intestine enriched?
"""),
code("""
activated_fh[activated_fh %in% liver_genes]
"""),
code("""
activated_fh[activated_fh %in% intestine_genes]
"""),
md("""
Is this more than expected by chance? We can test this using a **hypergeometric test**. Imagine an urn containing N balls (all genes), of which m are white (liver genes) and N - m are black (other genes). If we randomly draw k balls (the activated genes), how likely is it to draw q or more white balls?

- N: the total number of genes (our *universe*)
- m: the number of liver genes
- k: the number of activated genes
- q: the number of activated liver genes

The expected number of white balls is $k \\cdot \\frac{m}{N}$.

The choice of the universe is important. The HPA only contains protein coding genes, so we will use all protein coding genes in our count matrix (i.e. genes that are present in the HPA) as universe.
"""),
code("""
universe <- intersect(rownames(cts), hpa$Gene)
length(universe)
"""),
code("""
hypergeometric_test <- function(activated, tissue_genes, universe) {
    activated    <- intersect(activated, universe)
    tissue_genes <- intersect(tissue_genes, universe)
    N <- length(universe)
    m <- length(tissue_genes)
    k <- length(activated)
    q <- length(intersect(activated, tissue_genes))
    data.frame(
        observed = q,
        expected = k * m / N,
        # P(X >= q) = 1 - P(X <= q - 1)
        pvalue   = phyper(q - 1, m, N - m, k, lower.tail = FALSE))
}
"""),
code("""
enrichment_fh <- rbind(
    cbind(tissue = "liver",     hypergeometric_test(activated_fh, liver_genes, universe)),
    cbind(tissue = "intestine", hypergeometric_test(activated_fh, intestine_genes, universe)))
enrichment_fh
"""),
md("""
Let's plot this, similar to Figure 1B of the paper.
"""),
code("""
plot_enrichment <- function(enrichment, title) {
    d <- rbind(
        data.frame(tissue = enrichment$tissue, type = "expected", n = enrichment$expected),
        data.frame(tissue = enrichment$tissue, type = "observed", n = enrichment$observed))
    d$tissue <- factor(d$tissue, levels = c("liver", "intestine"))
    ggplot(d, aes(x = tissue, y = n, fill = type)) +
        geom_col(position = "dodge") +
        scale_fill_manual(values = c("expected" = "gray", "observed" = "darkred")) +
        ylab("number of activated genes") +
        ggtitle(title) +
        theme_bw()
}
options(repr.plot.width = 5, repr.plot.height = 4)
plot_enrichment(enrichment_fh, "FOXA1-HNF4A")
"""),
md("""
We find 88 activated liver genes, while we only expect ~ 14 genes by chance. Similarly, we find 30 intestinal genes, while we expect ~ 6 genes. Both enrichments are highly significant. This is very similar to Figure 1B of the paper (91 liver genes, 18 expected and 38 intestinal genes, 9 expected).

So, co-expression of FOXA1 and HNF4A in a blood cell line is sufficient to activate a large number of liver specific genes!
"""),
md("""
### FOXA1 and HNF4A separately

According to the pioneer factor hypothesis, FOXA1 (the pioneer factor) is needed to open chromatin, after which HNF4A (the non-pioneer factor) can bind and activate genes. If this is true, we expect that HNF4A alone activates (almost) no liver genes.

Let's repeat the analysis for the FOXA1 and HNF4A lines. We put all steps from above in a function.
"""),
code("""
run_deseq <- function(line) {
    samples <- rownames(coldata)[coldata$line == line]
    dds <- DESeqDataSetFromMatrix(countData = cts[, samples],
                                  colData = coldata[samples, ],
                                  design = ~ condition)
    dds$condition <- factor(dds$condition, levels = c("uninduced","induced"))
    keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
    dds <- dds[keep,]
    dds <- DESeq(dds, quiet = TRUE)
    resLFC <- lfcShrink(dds, coef="condition_induced_vs_uninduced", type = "apeglm", quiet = TRUE)
    list(dds = dds, resLFC = resLFC)
}
"""),
code("""
results_per_line <- list(
    FOXA1      = run_deseq("FOXA1"),
    HNF4A      = run_deseq("HNF4A"),
    FOXA1HNF4A = list(dds = dds, resLFC = resLFC))   # we already did this one
"""),
md("""
First, let's check the transgenes again.
"""),
code("""
do.call(rbind, lapply(names(results_per_line), function(line) {
    df <- as.data.frame(results_per_line[[line]]$resLFC)
    data.frame(line = line, gene = c("FOXA1", "HNF4A"),
               df[c("FOXA1", "HNF4A"), c("baseMean", "log2FoldChange", "padj")], row.names = NULL)
}))
"""),
md("""
- HNF4A is strongly induced in the HNF4A and FOXA1-HNF4A lines. In the FOXA1 line HNF4A is not expressed, the gene was even removed by our prefiltering step (hence the NA values).
- The counts of the endogenous FOXA1 gene go down upon induction in the FOXA1 and FOXA1-HNF4A lines (remember: the reads of the codon optimized transgene do not map to the FOXA1 gene).

**Question**: Can you think of an explanation why the expression of endogenous *FOXA1* goes down when the FOXA1 transgene is expressed?

Now let's get the activated genes for each line and test for liver and intestine enrichment.
"""),
code("""
activated <- lapply(results_per_line, function(x) get_activated_genes(x$dds, x$resLFC))
sapply(activated, length)
"""),
code("""
enrichment <- do.call(rbind, lapply(names(activated), function(line) {
    rbind(
        cbind(line = line, tissue = "liver",     hypergeometric_test(activated[[line]], liver_genes, universe)),
        cbind(line = line, tissue = "intestine", hypergeometric_test(activated[[line]], intestine_genes, universe)))
}))
enrichment
"""),
code("""
options(repr.plot.width = 12, repr.plot.height = 4)
d <- rbind(
    data.frame(line = enrichment$line, tissue = enrichment$tissue, type = "expected", n = enrichment$expected),
    data.frame(line = enrichment$line, tissue = enrichment$tissue, type = "observed", n = enrichment$observed))
d$tissue <- factor(d$tissue, levels = c("liver", "intestine"))
d$line   <- factor(d$line, levels = c("FOXA1", "HNF4A", "FOXA1HNF4A"))
ggplot(d, aes(x = tissue, y = n, fill = type)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("expected" = "gray", "observed" = "darkred")) +
    facet_wrap(~ line) +
    ylab("number of activated genes") +
    theme_bw()
"""),
md("""
- FOXA1 alone activates liver genes (26 observed, 10 expected), but no intestinal genes.
- HNF4A alone activates both liver (46 observed, 13 expected) and intestinal genes (28 observed, 6 expected).
- FOXA1-HNF4A activates the most liver genes.

This is in line with Figure 2A and 2B of the paper (FOXA1: 33 liver genes, intestine not significant, HNF4A: 47 liver genes, intestine significant).

**This is not what the pioneer factor hypothesis predicts!** HNF4A, the "non-pioneer" factor, can activate liver genes on its own, even more than FOXA1, the "pioneer" factor.
"""),
md("""
### Do FOXA1 and HNF4A activate the same liver genes?

Let's compare the activated liver genes between the three lines. An [UpSet plot](https://upset.app/) shows the size of each intersection: a dot means the gene is activated in that line, connected dots mean that the genes are activated in multiple lines.
"""),
code("""
library(UpSetR)
activated_liver <- lapply(activated, function(x) intersect(x, liver_genes))
sapply(activated_liver, length)
"""),
code("""
options(repr.plot.width = 7, repr.plot.height = 5)
upset(fromList(activated_liver), order.by = "freq", text.scale = 1.5)
"""),
md("""
- Only few liver genes are activated by both FOXA1 alone and HNF4A alone (6 + 3 genes): FOXA1 and HNF4A mostly activate *different* liver genes (Figure 2C-D of the paper).
- Most liver genes activated by FOXA1 or HNF4A alone are also activated by FOXA1-HNF4A.
- 42 liver genes are only activated when both factors are expressed together.

Genes that are activated by FOXA1-HNF4A but neither by FOXA1 nor HNF4A alone need *both* factors (Figure 4A of the paper):
"""),
code("""
need_both <- setdiff(activated_liver$FOXA1HNF4A, union(activated_liver$FOXA1, activated_liver$HNF4A))
length(need_both)
need_both
"""),
md("""
Let's visualize the expression of all activated liver genes in all samples using a heatmap. We use the variance stabilized counts of all samples (`vsd_all`, computed above).
"""),
code("""
options(repr.plot.width = 8, repr.plot.height = 12)
all_activated_liver <- unique(unlist(activated_liver))
all_activated_liver <- all_activated_liver[all_activated_liver %in% rownames(vsd_all)]
gene_annotation <- data.frame(
    row.names  = all_activated_liver,
    FOXA1      = ifelse(all_activated_liver %in% activated_liver$FOXA1, "yes", "no"),
    HNF4A      = ifelse(all_activated_liver %in% activated_liver$HNF4A, "yes", "no"),
    FOXA1HNF4A = ifelse(all_activated_liver %in% activated_liver$FOXA1HNF4A, "yes", "no"))
ann_colors <- list(FOXA1 = c(yes = "black", no = "white"), HNF4A = c(yes = "black", no = "white"), FOXA1HNF4A = c(yes = "black", no = "white"))
pheatmap(assay(vsd_all)[all_activated_liver, ],
         scale = "row", cluster_cols = FALSE,
         annotation_col = coldata[, c("line", "condition")],
         annotation_row = gene_annotation, annotation_colors = ann_colors,
         fontsize_row = 5)
options(repr.plot.width = 7, repr.plot.height = 5)
"""),
md("""
The genes cluster in groups:

- genes activated by HNF4A (and FOXA1-HNF4A), for example the apolipoproteins APOA2 and APOC3 and SERPINA1.
- genes activated by FOXA1 (and FOXA1-HNF4A), for example ALB and the fibrinogens FGB and FGG.
- genes activated by both FOXA1 and HNF4A alone, for example APOA1 and ARG1.
- genes that are only activated in the FOXA1-HNF4A line, for example FGA, PLG (plasminogen) and F2 (prothrombin). These genes need both factors.
"""),
md("""
### What is the function of the activated liver genes? Gene ontology (GO) enrichment

[Gene Ontology](http://geneontology.org/) (GO) is a database that describes the function of genes using a controlled vocabulary of terms. These terms are organized in three categories:

- Biological Process (BP), e.g. "cholesterol metabolic process"
- Molecular Function (MF), e.g. "DNA binding"
- Cellular Component (CC), e.g. "nucleus"

Similar to the liver enrichment above, we can test for each GO term whether it is enriched in a list of genes compared to a universe of genes using a hypergeometric test. This is implemented in the function `enrichGO` of the package [clusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/). Because we test thousands of GO terms, we need to correct for multiple testing.

In the paper (Supplementary files 1 and 2), the liver genes activated by FOXA1 were enriched for fibrinolysis and complement activation, while the liver genes activated by HNF4A were enriched for cholesterol import and lipoprotein remodeling. Let's see if we find the same.

As universe we use all genes that were tested in the corresponding DESeq2 analysis (after filtering).
"""),
code("""
suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
})
"""),
code("""
run_go <- function(genes, universe) {
    enrichGO(gene          = genes,
             universe      = universe,
             OrgDb         = org.Hs.eg.db,
             keyType       = "SYMBOL",   # our genes are gene symbols (gene names)
             ont           = "BP",       # biological process
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.2,
             minGSSize     = 5)
}
"""),
code("""
go_foxa1 <- run_go(activated_liver$FOXA1, rownames(results_per_line$FOXA1$dds))
head(as.data.frame(go_foxa1)[, c("Description", "GeneRatio", "BgRatio", "p.adjust", "geneID")], 15)
"""),
code("""
dotplot(go_foxa1, showCategory = 15, font.size = 10) + theme(plot.margin = margin(5, 5, 5, 60)) + ggtitle("FOXA1 activated liver genes")
"""),
code("""
go_hnf4a <- run_go(activated_liver$HNF4A, rownames(results_per_line$HNF4A$dds))
head(as.data.frame(go_hnf4a)[, c("Description", "GeneRatio", "BgRatio", "p.adjust", "geneID")], 15)
"""),
code("""
dotplot(go_hnf4a, showCategory = 15, font.size = 10) + theme(plot.margin = margin(5, 5, 5, 60)) + ggtitle("HNF4A activated liver genes")
"""),
md("""
As in the paper:

- the liver genes activated by FOXA1 are enriched for terms related to blood coagulation, **fibrinolysis** (the breakdown of blood clots, e.g. KNG1, FGB, FGG, SERPINF2) and **complement activation** (part of the innate immune system).
- the liver genes activated by HNF4A are enriched for terms related to lipid metabolism, **(plasma) lipoprotein particle remodeling** (e.g. the apolipoproteins) and bile acid metabolism, but also complement activation.

So FOXA1 and HNF4A each activate a different part of the liver gene expression program.
"""),
md("""
We can also look at all activated genes (not only the liver genes) of the FOXA1-HNF4A line.
"""),
code("""
go_fh <- run_go(activated$FOXA1HNF4A, rownames(dds))
options(repr.plot.width = 10, repr.plot.height = 7)
dotplot(go_fh, showCategory = 20, font.size = 10) + theme(plot.margin = margin(5, 5, 5, 60)) + ggtitle("FOXA1-HNF4A activated genes")
"""),
md("""
The activated genes are enriched for typical liver (and intestine) functions: transport of organic acids, steroid metabolism, response to xenobiotics (detoxification of drugs and toxins), blood coagulation, digestion and intestinal lipid absorption.
"""),
md("""
### Comparison with the paper

| | this notebook | paper |
|---|---|---|
| FOXA1-HNF4A: activated liver genes (expected) | 88 (14) | 91 (18) |
| FOXA1-HNF4A: activated intestine genes (expected) | 30 (6) | 38 (9) |
| FOXA1: activated liver genes | 26 | 33 |
| FOXA1: intestine enrichment | not significant | not significant |
| HNF4A: activated liver genes | 46 | 47 |
| HNF4A: intestine enrichment | significant | significant |
| liver genes only activated by FOXA1-HNF4A | 42 | 31 |
| GO FOXA1 liver genes | fibrinolysis, complement activation | fibrinolysis, complement activation |
| GO HNF4A liver genes | lipoprotein particle remodeling | cholesterol import, lipoprotein remodeling |

Our numbers are not exactly the same as in the paper. This is expected, because:

- we used a different pipeline (STAR + featureCounts instead of Salmon), and a different version of the genome annotation.
- we required a padj < 0.05, the paper only used a log2 fold change threshold.
- we used shrunken log2 fold changes.
- we used a newer version of the Human Protein Atlas.
- the universe used for the hypergeometric test is not described in the paper.

Nevertheless, the conclusions are the same!
"""),
md("""
## Conclusion

1. Expression of FOXA1 and HNF4A in K562 cells (a blood cell line) activates more than 2000 genes, among which many liver and intestine specific genes.
2. Both FOXA1 alone and HNF4A alone can activate liver genes. HNF4A alone also activates intestinal genes.
3. FOXA1 and HNF4A activate largely different sets of liver genes, and a subset of liver genes requires both factors.

According to the pioneer factor hypothesis only the pioneer factor FOXA1 should be able to activate silent genes on its own. The RNA-seq data shows that the "non-pioneer" factor HNF4A can do so as well. Using ATAC-seq and CUT&Tag data (not analysed here), the authors show that both factors can bind to and open closed chromatin. They conclude that "pioneer activity" is not a property of a special class of transcription factors, but depends on the binding affinity of the transcription factor and the genomic context.
"""),
md("""
## Exercises

1. Change the threshold for "activated" genes (e.g. `max_uninduced = 10` or `log2FC_thr = 2`). How does this change the number of activated liver genes and the enrichment?
2. Repeat the hypergeometric test for other tissues in the HPA (e.g. `"brain"`, `"bone marrow"`, `"pancreas"`). Are these tissues also enriched? What do you expect?
3. Which genes are *downregulated* upon FOXA1-HNF4A induction? Run a GO enrichment on these genes. Do you see a decrease in blood cell related genes?
4. Instead of analysing each line separately, make a single DESeqDataSet with all 18 samples and design `~ group` and use `results(dds_all, contrast = c("group", "FOXA1_induced", "FOXA1_uninduced"))`. Do you get the same results? (Hint: what is the difference in how the dispersion is estimated?)
5. Plot the normalized counts of the genes `ALB`, `APOA2` and `AMDHD1` across all 18 samples. Which factor(s) activate these genes?
"""),
]

if __name__ == "__main__":
    nb = new_notebook(cells, "ir-foxa1_hnf4a", "R (FoxA1/Hnf4a RNA-seq)", "R")
    if "--execute" in sys.argv:
        nb = execute(nb, cwd=TOPIC)
    elif "--markdown-only" in sys.argv:
        nb = sync_markdown(OUT, nb)
    nbf.write(nb, OUT)
    print("written", OUT)
