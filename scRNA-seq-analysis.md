scRNA-seq analysis of muscle sample
================
Paola Maragno, Alberto Pettenella
17-06-2022

``` r
library(dplyr)
library(Seurat)
library(patchwork)
library("stringr")
```

# Data presentation

In this project we will analyze single cell RNA-seq data coming from
2737 cells extracted from a muscle sample of a 3 month mouse (strain:
C57BL/6 NIA).

**The approach that was used to perform scRNA-seq is 10x Genomics
Chromium which is based on the GemCode technology:** first the cells are
combined with reagents in one channel of a microfluidic chip and gel
beads from another channel to form GEMs (Gel bead in EMulsion). Each gel
bead is functionalized with barcoded oligonucleotides which consists of:
adapters and primers, a barcode, a unique molecular identifier (UMI) and
an anchored 30 bp oligo-dT to prime polyadenylated RNA transcripts.

After encapsulation of cells in GEMs, the gel beads dissolve and release
their oligonucleotides for reverse transcription of polyadenylated RNAs:
each resulting cDNA molecule contains its own UMI and a shared barcode
for all the cDNAs coming from the same GEM, and ends with a template
switching oligo at the 3′ end. Next, each GEM is broken and barcoded
cDNAs are pooled for PCR amplification.

Eventually, amplified cDNAs are fragmented and undergo NG sequencing.
Read1 tags the cDNA sequence at the 3’-polyA end, while Read2 captures
the UMI and the barcode. Reads with the same UMI are replicates of the
same initial cDNA fragment, while reads with the same barcode and
different UMI are replicates of different cDNAs coming from the same
cell.

All the sequences of the reads obtained through this process performed
on the cells of the muscle sample are downloaded from PanglaoDB.
**According to PanglaoDB the cell types present in this sample are:**

| **Cell type**     | **Number of cells** |
|-------------------|---------------------|
| B cells           | 160                 |
| Endothelial cells | 842                 |
| Fibroblast        | 739                 |
| Macrophages       | 135                 |
| Myoblasts         | 26                  |
| Pericytes         | 82                  |
| Schwann cell      | 16                  |
| T memory cells    | 128                 |
| Unknown           | 14                  |

While **according to the article from which the data were collected
(<https://pubmed.ncbi.nlm.nih.gov/30283141/>) the cell types present in
this sample are:**

-   B cells
-   T cells
-   Endothelial cells
-   Macrophages
-   Mesenchymal stem cells (that corresponds to the Myoblasts)
-   Skeletal muscle satellite cells

# Data loading and pre-processing

First we upload the **Sm object that is a count table** containing as
columns the barcode of each single cell that has been isolated in a
droplet, while as rows the names of the genes which expression was
measured in each cell.

``` r
load('SRA653146_SRS3044259.sparse.RData')
```

First we check the names of the genes in the count table

``` r
head(rownames(sm))
```

    ## [1] "00R_AC107638.2_ENSMUSG00000111425.1" "0610009O20Rik_ENSMUSG00000024442.5" 
    ## [3] "1010001N08Rik_ENSMUSG00000097222.1"  "1110020A21Rik_ENSMUSG00000097047.1" 
    ## [5] "1700007L15Rik_ENSMUSG00000097318.1"  "1700015F17Rik_ENSMUSG00000079666.8"

We can see that the name of each gene is made of the gene symbol plus
the Ensembl ID separated by ’\_’. Consequently the following code
removes the unnecessary part and returns only the gene symbols:

``` r
rownames(sm) <- lapply(rownames(sm), FUN = function(x) {
      if (str_count(x, "_") != 1) {
        a <- strsplit(x, '_')
        x <- paste(a[[1]][1],a[[1]][2],sep='_')
  }
  else {
        x <-strsplit(x, '_')[[1]][1]}
  })
  
head(rownames(sm))
```

    ## [1] "00R_AC107638.2" "0610009O20Rik"  "1010001N08Rik"  "1110020A21Rik" 
    ## [5] "1700007L15Rik"  "1700015F17Rik"

Now we can **initialize the Seurat object with the raw data** contained
in the just modified count table so that the row names will be the
correct ones. The parameters that are used when creating this object
allow to filter out from the count table all the cells in which there
are less than 200 genes expressed and also all the genes that are
expressed in less than 3 cells

``` r
pbmc <- CreateSeuratObject(
     sm,
     project = "scRNA_seq_muscle",
     min.cells = 3,
     min.features = 200
)

head(rownames(pbmc))
```

    ## [1] "00R-AC107638.2" "0610009O20Rik"  "1010001N08Rik"  "1110020A21Rik" 
    ## [5] "1700007L15Rik"  "1700022N22Rik"

We can see that the row names are just the gene symbols.

``` r
pbmc
```

    ## An object of class Seurat 
    ## 20372 features across 2737 samples within 1 assay 
    ## Active assay: RNA (20372 features, 0 variable features)

2737 cells will be analyzed (the samples) using the expression value of
20372 genes (the features) coming from them.

As said the names of the columns of the count table are the barcodes of
the single cells that were analyzed

``` r
head(colnames(sm))
```

    ## [1] "AAACCTGAGCAGATCG" "AAACCTGAGTACTTGC" "AAACCTGGTCTAGTGT" "AAACCTGTCCTGCTTG"
    ## [5] "AAACCTGTCGCCTGAG" "AAACGGGAGAATTCCC"

Let’s examine the count matrix for three different genes in the first 30
cells

``` r
sm[c("Fabp4", "Gsn", "Fmod"), 1:30]
```

    ## 3 x 30 sparse Matrix of class "dgCMatrix"
    ##                                                                                
    ## Fabp4 177 .   . . 16 1 . 22 14 38 2 . 30 1   2 5 . 1 1 19 29 . .  .   . . 2 . 7
    ## Gsn     1 1 123 .  . 1 3  .  .  1 . 1  3 . 333 2 . . 1  1  1 2 .  9 199 1 . 2 1
    ## Fmod    . .   . .  . . .  .  .  . . .  1 .   . . . . 1  .  . 1 . 58   1 . . 1 .
    ##         
    ## Fabp4 66
    ## Gsn    1
    ## Fmod   .

‘.’ indicates that the gene is not expressed in a specific cell

# Cell quality control

The **main quality control parameters** that we have to control are:

-   the number of unique genes detected in each cell: if too few the
    cells may be of low-quality or empty droplets; if too many it can be
    consequence of cell doublets or multiplets;
-   the percentage of reads that map to the mitochondrial genome
    corresponding to “suffering” or “dead” cells, damaged during
    separation/processing;
-   not too many reads mapping on ribosomal protein genes since they
    will “eat up” a lot of reads because highly expressed.

To evaluate these aspects we add inside the Seurat object the
information about the percentage of counts originating from
mitochondrial genes or ribosomal protein genes:

``` r
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
pbmc[["percent.rbp"]] <- PercentageFeatureSet(pbmc, pattern = "^Rp[Ls]")
```

``` r
head(pbmc@meta.data, 5)
```

    ##                        orig.ident nCount_RNA nFeature_RNA percent.mt
    ## AAACCTGAGCAGATCG scRNA_seq_muscle       8901         2576  1.0448264
    ## AAACCTGGTCTAGTGT scRNA_seq_muscle       5203         1823  0.9225447
    ## AAACCTGTCGCCTGAG scRNA_seq_muscle       1284          706  1.6355140
    ## AAACGGGCACACCGAC scRNA_seq_muscle       4111         1317  0.9486743
    ## AAACGGGCATCTGGTA scRNA_seq_muscle       1203          734  2.0781380
    ##                  percent.rbp
    ## AAACCTGAGCAGATCG    8.954050
    ## AAACCTGGTCTAGTGT   10.724582
    ## AAACCTGTCGCCTGAG    5.218069
    ## AAACGGGCACACCGAC   21.722209
    ## AAACGGGCATCTGGTA    2.410640

**Meta data object contained in Seurat object collects useful
information for each cell** (that is identified by its barcode): the
overall number of reads coming from each cell, the number of genes
(features) that were found to be transcribed with at least one read in
that particular cell, the percentage of reads mapping on mitochondrial
genes and the percentage of those codifying for ribosomal proteins.

Visualize QC metrics as violin plots

``` r
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Each dot is a cell; if we want to see the same distributions without
dots:

``` r
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4, pt.size=0)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

We can see that the number of genes found transcribed in each cell is
distributed between 0 and 4000 and then there is the characteristic
single line that represents the outliers. The number of reads coming
from each cell varies between 0 and 10500 and then there are outliers
with higher values; the number of reads mapping on mitochondrial genes
is concentrated around low percentages (lower than 5%) with lot of
outliers. Eventually we can notice that all the cells have a percentage
of reads mapping on rbp genes around 0-30% without outliers: this is
very good!

We can check if the different parameters are correlated with one another

``` r
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

We visualize the **distribution of reads mapping on mitochondrial genes
with respect to the overall number of reads in each cell** and we can
assess that there is no correlation between the two, meaning that having
lot of reads doesn’t mean having lot/few reads mapping on mitochondrial
genes.

``` r
plot2
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

We can see that there is a strong positive **correlation between how
many reads derive from a cell and the number of genes found to be
transcribed in that cell.** It is important that the number of genes
coming from a cell is not too high otherwise there is the risk of being
analyzing a doublet or multiplet.

``` r
plot3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.rbp")
plot3
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

From the plot is visible **the lack of correlation between the number of
reads mapping on rbp genes and the total number of reads coming from a
cell:** this means that having a lot of reads doesn’t mean that lot of
them necessarily map on ribosomal protein genes.

To sum up what we can say from these visualizations is that **there is a
correlation only between the number of reads and the number of genes
detected in each cell.** On the basis of the first violin plot, we
decide to **keep for the further analysis only the cells with a number
of expressed genes between 200 and 4000 - so we remove the outliers that
express too many genes since they are probably doublets - and we remove
also all those cells which are outliers for the percentage of reads
mapping on mitochondrial genes - since they may be damaged cells**

``` r
pbmc <- subset(pbmc, subset = (nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 4))
```

Let us see how many samples remain:

``` r
pbmc
```

    ## An object of class Seurat 
    ## 20372 features across 2602 samples within 1 assay 
    ## Active assay: RNA (20372 features, 0 variable features)

135 cells have been throw away.

# Normalizing the data

For 10x platform the advice is to not normalize the data: this because
our aim is to cluster cells of the same cell type and to do so what is
more important to consider is the “signature” of each cell (so if that
cell expresses or not a gene) rather than the level at which the cell
expresses it. In other words, in the subsequent analyses what makes the
difference is how a gene changes its expression across the cells, rather
than its actual expression values.

For this reason for **10x data it is advised to normalize the data just
by multiplying the number of counts of that gene in that cell by 10000
to make it more human readable and then compute the logarithm of the
normalized counts.** After normalization is also advised to **“scale”
the log2-counts per million of each gene across the cells** to shift the
expression of each gene so that the mean expression across cells is 0
and its variance 1: this step gives equal weight in downstream analyses,
so that highly-expressed genes do not dominate.

``` r
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```

Now, in the Seurat object we can find the original counts using the
following code: <pbmc@assays>$<RNA@counts>

And also the normalized ones: <pbmc@assays>$<RNA@data>

Seurat contains a **pre-computed list of cell cycle specific genes**: so
genes known to be expressed in particular stages of cell cycle.
**CellCycleScoring function, according to which of the marker genes of S
and G2M cycle phases each cell expresses, guesses in which phase of the
cell cycle each cell is**

``` r
CellCycleScoring(pbmc, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> pbmc
```

``` r
head(pbmc)
```

    ##                        orig.ident nCount_RNA nFeature_RNA percent.mt
    ## AAACCTGAGCAGATCG scRNA_seq_muscle       8901         2576  1.0448264
    ## AAACCTGGTCTAGTGT scRNA_seq_muscle       5203         1823  0.9225447
    ## AAACCTGTCGCCTGAG scRNA_seq_muscle       1284          706  1.6355140
    ## AAACGGGCACACCGAC scRNA_seq_muscle       4111         1317  0.9486743
    ## AAACGGGCATCTGGTA scRNA_seq_muscle       1203          734  2.0781380
    ## AAACGGGCATGTCTCC scRNA_seq_muscle        679          473  1.6200295
    ## AAACGGGGTGTCTGAT scRNA_seq_muscle       1339          686  1.7176998
    ## AAAGATGAGCTCAACT scRNA_seq_muscle       3463         1513  0.6064106
    ## AAAGATGAGTGATCGG scRNA_seq_muscle       8608         2436  0.8480483
    ## AAAGATGCACTGTCGG scRNA_seq_muscle       3419         1390  0.5264697
    ##                  percent.rbp      S.Score    G2M.Score Phase        old.ident
    ## AAACCTGAGCAGATCG    8.954050 -0.020789460 -0.045433394    G1 scRNA_seq_muscle
    ## AAACCTGGTCTAGTGT   10.724582 -0.025739636 -0.052994242    G1 scRNA_seq_muscle
    ## AAACCTGTCGCCTGAG    5.218069  0.025923669  0.000264755     S scRNA_seq_muscle
    ## AAACGGGCACACCGAC   21.722209 -0.051401347  0.053023077   G2M scRNA_seq_muscle
    ## AAACGGGCATCTGGTA    2.410640 -0.020829601  0.079728154   G2M scRNA_seq_muscle
    ## AAACGGGCATGTCTCC    6.480118  0.036784352  0.003273698     S scRNA_seq_muscle
    ## AAACGGGGTGTCTGAT   10.978342  0.043866795 -0.055943038     S scRNA_seq_muscle
    ## AAAGATGAGCTCAACT    3.869477  0.005693484 -0.058024215     S scRNA_seq_muscle
    ## AAAGATGAGTGATCGG    5.425186 -0.045062137 -0.070829671    G1 scRNA_seq_muscle
    ## AAAGATGCACTGTCGG    3.305060 -0.023850186 -0.024590313    G1 scRNA_seq_muscle

Now the phase of the cell cycle in which the cell is estimated to be is
added to the Seurat object: G1 means the cell is not cycling, S means it
is in the replicative phase, while G2M means it is in the transition
between phase G2 and M.

**Each cell is a point in a n-dimensional space**, where n is the number
of genes for which we have the expression value in each cell. The closer
two points are, the more similar are the transcriptomes of the
corresponding cells and so these points should be clustered together
since they belong to the same cell type.

Since the dimensions are too many for further processing and the
expression of lot of genes in each cell is zero, it is better to keep
only the subset of genes with the greatest variability of expression
across cells. **The default method -vst- estimates the mean-variance
relationship of each gene and chooses the 2000 genes with the highest
variance**

``` r
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
```

In this way we have reduced the 20372-dimensional space in a
2000-dimensional one: all the genes with 0 expression in all the cells,
all those with “random” variations (genes peaking in just a few, e.g. \<
10, cells) and all housekeeping genes with signature showing them to be
expressed by most of the cells were thrown away.

Identify the 10 most highly variable genes

``` r
top10 <- head(VariableFeatures(pbmc), 10)
```

Plot the 2000 variable features with and without labels of the top 10
most variable ones

``` r
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
plot2
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

The plots show in red the 2000 most variable features - on the bottom
there is also the label of the 10 most variable ones - while in black
there are the 18372 not variable genes that won’t be taken in account in
the further analysis.

Now **the scaling of the genes can be performed**

``` r
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

Inside the Seurat object we have now also the normalized and scaled
counts: for example here we can see the normalized and scaled counts of
Fabp4 gene in the first 6 cells

``` r
head(pbmc[["RNA"]]@scale.data["Fabp4",])
```

    ## AAACCTGAGCAGATCG AAACCTGGTCTAGTGT AAACCTGTCGCCTGAG AAACGGGCACACCGAC 
    ##        1.4178103       -0.9539746        1.2098932       -0.9539746 
    ## AAACGGGCATCTGGTA AAACGGGCATGTCTCC 
    ##        1.3805048        1.4339401

What is important to check is if the similarity between cells is due to
the fact that they are in the same cell cycle phase: in some cases this
is something we want - for example in case the sample is made of non
proliferating and proliferating cells so the fact that the cells are
cycling is something that makes the difference. Since **in this case the
information related to the phase of the cell cycle in which each cell is
is something not relevant to detect the similarity between cells - but
we want to evaluate it only according to their transcriptom profile - we
now must check if it is necessary to remove the cell cycle effect.**

The best way to check it is visualizing the cells distribution in the
space using a technique of dimensionality reduction: the PCA. Notice
that **PCA is performed only on the most variable features**

``` r
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```

    ## PC_ 1 
    ## Positive:  Col6a1, Dcn, Col1a2, Serping1, Col1a1, Gsn, Igfbp6, Col3a1, Col6a2, Pcolce 
    ##     Serpinf1, Mmp2, Fstl1, Clec3b, Timp2, Rnase4, Ogn, Mfap5, Nbl1, Il11ra1 
    ##     Axl, Col5a2, Dpt, Islr, Htra3, Col6a3, Bgn, Lum, Fndc1, Tnxb 
    ## Negative:  Tmsb4x, Fabp4, Cd52, Coro1a, Laptm5, Ptpn18, Rac2, Gpihbp1, Aqp1, Cd53 
    ##     Gmfg, H2-Aa, Cd37, Tinagl1, Egfl7, Arhgap45, Cd74, H2-Eb1, Ctss, H2-Ab1 
    ##     C1qtnf9, RP23-310J6.1, Arhgdib, Cytip, H2-DMa, Ctla2a, Ppp1r16b, Id1, Apold1, H2-DMb2 
    ## PC_ 2 
    ## Positive:  Fabp4, Hspb1, Tinagl1, Gpihbp1, Aqp1, Id1, Apold1, RP23-310J6.1, Tm4sf1, Egfl7 
    ##     Rgcc, C1qtnf9, Sox17, Cxcl12, Tcf15, Timp4, Rnd1, Cldn5, Iigp1, Ctla2a 
    ##     Fabp5, Grrp1, Icam2, Cdc42ep3, Tsc22d1, Isg15, RP23-303F24.3, Thbd, Hes1, RP24-171J16.5 
    ## Negative:  Ctss, Fcer1g, Tyrobp, Alox5ap, Lyz2, Laptm5, Cd74, H2-Eb1, Ccl9, H2-DMa 
    ##     H2-Ab1, H2-Aa, Spi1, Fcgr3, Cd53, Cd52, Cd68, Ccl6, H2-DMb1, C1qc 
    ##     C1qb, Ly86, Ptafr, Csf1r, C1qa, Cd83, Pld4, Wfdc17, Ctsc, Coro1a 
    ## PC_ 3 
    ## Positive:  Cst3, Fabp4, Ier3, C1qa, C1qb, C1qc, Csf1r, Tm4sf1, Fcgr3, Aqp1 
    ##     Lyz2, Tinagl1, Fcer1g, Apold1, Dusp1, Pf4, Ccl9, Hspb1, Mrc1, Alox5ap 
    ##     Crip1, Atf3, Tmsb4x, Cd68, Rnd1, Ccl6, Id1, C5ar1, Cd14, C3ar1 
    ## Negative:  Chchd10, Satb1, Ptprcap, Pgam2, Ckm, RP24-146B4.2, RP24-426K19.3, Cd79a, Atp2a1, Cox6a2 
    ##     Myl1, Ltb, Tcap, Tnni2, Actn3, Acta1, Pvalb, Eno3, Myoz1, Car3 
    ##     Tnnt3, Ccr7, Cox8b, Cd79b, Ighm, RP23-32C18.6, Ighd, Mylpf, Limd2, Rac2 
    ## PC_ 4 
    ## Positive:  Dpt, Cd34, Clec3b, Ifi205, Dpep1, Cd248, Scara5, Ccl11, Cd55, Pi16 
    ##     Cadm3, Mfap5, Efemp1, Plac9c, Efhd1, Ace, Entpd2, Pdgfra, Ptx3, Pla1a 
    ##     Mnda, Ifi204, Tmsb4x, Col14a1, Pcsk6, Emilin2, Nid1, Lrrn4cl, Ms4a4d, Plac8 
    ## Negative:  Cilp2, Tnmd, Fmod, Chad, Thbs4, Col11a1, Comp, Ptx4, Col12a1, Cpxm2 
    ##     Kera, Atp2a1, Ckm, RP24-426K19.3, Cox6a2, Myl1, Pgam2, Tcap, Eno3, Itgbl1 
    ##     Tpm2, Tpm1, Scx, Cox8b, Acta1, Col8a2, Mylpf, Myoz1, Mfap4, AC159092.1 
    ## PC_ 5 
    ## Positive:  Cilp2, Tnmd, Fmod, Thbs4, Comp, Chad, Cpxm2, Col11a1, Col12a1, Ptx4 
    ##     Kera, Itgbl1, Mfap4, AC159092.1, Cilp, Col8a2, Ptgis, Col11a2, Angptl7, Wif1 
    ##     Tnc, Scx, Mkx, Clec11a, Lox, Olfml3, Emb, Kctd1, Wisp1, Timp1 
    ## Negative:  Ckm, RP24-426K19.3, Pgam2, Myl1, Atp2a1, Cox6a2, Tcap, Tnni2, Myoz1, Acta1 
    ##     Pvalb, Tnnt3, Mylpf, Eno3, Car3, Actn3, Cox8b, Smpx, Eef1a2, Adssl1 
    ##     Rpl3l, Cox7a1, Apobec2, Myot, Hspb6, Tmod4, Mybpc2, Tpm2, Ak1, Sh3bgr

We can see which are the genes that make the higher variability along
the different principal components. **The most variable features along
the first PCs are those that make the difference in discovering which
are the different cell types inside the sample**

``` r
print(pbmc[["pca"]], dims = 1:5, , nfeatures = 5)
```

    ## PC_ 1 
    ## Positive:  Col6a1, Dcn, Col1a2, Serping1, Col1a1 
    ## Negative:  Tmsb4x, Fabp4, Cd52, Coro1a, Laptm5 
    ## PC_ 2 
    ## Positive:  Fabp4, Hspb1, Tinagl1, Gpihbp1, Aqp1 
    ## Negative:  Ctss, Fcer1g, Tyrobp, Alox5ap, Lyz2 
    ## PC_ 3 
    ## Positive:  Cst3, Fabp4, Ier3, C1qa, C1qb 
    ## Negative:  Chchd10, Satb1, Ptprcap, Pgam2, Ckm 
    ## PC_ 4 
    ## Positive:  Dpt, Cd34, Clec3b, Ifi205, Dpep1 
    ## Negative:  Cilp2, Tnmd, Fmod, Chad, Thbs4 
    ## PC_ 5 
    ## Positive:  Cilp2, Tnmd, Fmod, Thbs4, Comp 
    ## Negative:  Ckm, RP24-426K19.3, Pgam2, Myl1, Atp2a1

**Col6a1** is found in macrophages, B cells, T cells, endothelial cells,
skeletal muscle satellite cells and mesenchymal stem cells; while
**Fabp4** is found in endothelial cells and macrophages. This is a
suggestion of the fact that one or more of these cell types are present
in the sample.

Plot the most variable genes on the first two PCs

``` r
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

And the projection of the cells in the first two principal components

``` r
DimPlot(pbmc, reduction = "pca")
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

Cells are plotted in the bidimensional space colored in a different way
depending on the cell cycle phase in which they are. If cell cycle is
the main source of difference between our cells we would see a cluster
of green/red/blue cells, meaning that the first 2 PC, that are those
that explain the greatest part of variability among the cells, partition
the cells according to the cell cycle and thus it is an important factor
in evaluating the cell similarity.

**In our plot the cells of different CC phases are mixed suggesting that
the cell cycle factor doesn’t need to be removed.**

# Clustering

To proceed in the further analysis we have to further reduce the
dimensionality of our data since it is not possible to process 2000
genes for each cell but we need to capture only those very meaningful
genes among them. **We use the PCA technique to perform dimensionality
reduction**: each PC combines information across a correlated feature
set. So the top principal components are those able to explain the
greatest variability in our data.

To evaluate which number of principal components is able to explain our
data in a satisfactory way we use the **Elbow plot** in which, for each
number of selected PCs, the overall SD of the data that remains not
explained is reported

``` r
ElbowPlot(pbmc, ndims=50)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

The rule is choosing the first number of PCs for which the overall SD of
the data is lower than 2 and where the plateau starts. From the plot the
first number of PC for which these two events happen is more or less 25.

There are different ways to perform clustering: the best one on 10X data
are **graph-based methods like the Louvain algorithm that returns the
clusters of data that are the best from a biological point of view.**
Graph-based clustering works in this way: cells are represented by nodes
in the graph and each cell is connected by an edge to its k-Nearest
Neighbors so the k cells that are more similar to it - so that are
closer to it in the space. Furthermore the edges between nodes are
“weighted” since each edge is weighted for the similarity value between
the cells that are connected by it. The similarity value is a mix of how
cells are far away in the Euclidean space plus how many similar cells
these two cells have in common.

**The best partition of the graph is the one in which the sum of the
similarity values between cells of the same cluster is the maximum,
while the sum of the similarity values between cells belonging to
different clusters should be the minimum. So we want to have the highest
number of edges connecting cells of the same cluster and the lowest,
ideally 0, number of edges that connect cells of different clusters.**

To compute in automatic way the best number of PCs we can use the
following code that keeps all the PC until all in all we are able to
explain the 75% of the overall standard deviation

``` r
pc.touse <- (pbmc$pca@stdev)^2
pc.touse <- pc.touse/sum(pc.touse)
pc.touse <- cumsum(pc.touse)[1:50]
pc.touse <- min(which(pc.touse>=0.75))
pc.touse
```

    ## [1] 16

Even if from the Elbow plot it seems that the best number of PCs is
around 25/30, from this code we can see that already using 16 PCs we are
able to explain great part of the variability inside the data. Since
using more PCs has the drawback of introducing more noise and so even
the result of the clustering will be noisier, it is better to consider a
lower number of PCs knowing that they are already able to deal with a
lot of variability.

With the function **FindNeighbors for each cell the k neighbors cells
are found:** these neighbors are the cells which transcriptom profile,
considering only the genes that make the higher variability along the
first 16 principal components, is more similar to the one of the cell
taken into consideration.

``` r
pbmc16_02 <- FindNeighbors(pbmc, dims = 1:16)
```

Once the neighbors of each cell are found the **FindClusters() function
implements the Louvain algorithm: we can specify a resolution parameter
that sets the ‘granularity’ of the downstream clustering.** Increasing
the resolution value we obtain a greater number of clusters. Seurat
authors found that setting this parameter between 0.4-1.2 typically
returns good results for single-cell datasets of around 3K cells: since
we have 2602 cells we also try values lower than 0.4.

``` r
pbmc16_02 <- FindClusters(pbmc16_02, resolution = 0.2)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2602
    ## Number of edges: 80813
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9586
    ## Number of communities: 12
    ## Elapsed time: 0 seconds

The last two commands, so the detection of the k neighbors of each cell
considering a space with n-dimensions, where n is the number of
principal components we select, and the identification of clusters with
the assignment of each cell to a cluster are repeated different times
trying different values of the two parameters - number of principal
components and resolution - till we achieve a satisfactory solution.

Once the clustering was performed, **in the pbmc object we can find the
information of the cluster to which each cell is associated**

``` r
head(pbmc16_02[[]],5)
```

    ##                        orig.ident nCount_RNA nFeature_RNA percent.mt
    ## AAACCTGAGCAGATCG scRNA_seq_muscle       8901         2576  1.0448264
    ## AAACCTGGTCTAGTGT scRNA_seq_muscle       5203         1823  0.9225447
    ## AAACCTGTCGCCTGAG scRNA_seq_muscle       1284          706  1.6355140
    ## AAACGGGCACACCGAC scRNA_seq_muscle       4111         1317  0.9486743
    ## AAACGGGCATCTGGTA scRNA_seq_muscle       1203          734  2.0781380
    ##                  percent.rbp     S.Score    G2M.Score Phase        old.ident
    ## AAACCTGAGCAGATCG    8.954050 -0.02078946 -0.045433394    G1 scRNA_seq_muscle
    ## AAACCTGGTCTAGTGT   10.724582 -0.02573964 -0.052994242    G1 scRNA_seq_muscle
    ## AAACCTGTCGCCTGAG    5.218069  0.02592367  0.000264755     S scRNA_seq_muscle
    ## AAACGGGCACACCGAC   21.722209 -0.05140135  0.053023077   G2M scRNA_seq_muscle
    ## AAACGGGCATCTGGTA    2.410640 -0.02082960  0.079728154   G2M scRNA_seq_muscle
    ##                  RNA_snn_res.0.2 seurat_clusters
    ## AAACCTGAGCAGATCG               8               8
    ## AAACCTGGTCTAGTGT               1               1
    ## AAACCTGTCGCCTGAG               0               0
    ## AAACGGGCACACCGAC               6               6
    ## AAACGGGCATCTGGTA               0               0

To look at the cluster IDs of the first 5 cells

``` r
head(Idents(pbmc16_02), 5)
```

    ## AAACCTGAGCAGATCG AAACCTGGTCTAGTGT AAACCTGTCGCCTGAG AAACGGGCACACCGAC 
    ##                8                1                0                6 
    ## AAACGGGCATCTGGTA 
    ##                0 
    ## Levels: 0 1 2 3 4 5 6 7 8 9 10 11

We **plot the cells in the space of the first two PCA components**

``` r
DimPlot(pbmc16_02, reduction = "pca")
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

The points, so the cells, in the PCA space are colored by cluster
membership. We can see that in some clusters the cells are very
concentrated but in other clusters they are spread around and this
because **some clusters are well explained by the first two PCs, while
others are better explained by other PCs.**

Another way to visualize the cells divided in clusters is
**“t-stochastic neighbor embedding”**: it attempts to find a
low-dimensional representation of the data that preserves the distances
between each point and its neighbors in the high-dimensional space -
that is the space made of the number of principal components that were
used to find the clusters. This method tends to display more “beautiful”
clusters in the data and allows to overcome the fact that the first two
PCs are not able to represent in the same way all the clusters.

``` r
pbmc16_02 <- RunTSNE(pbmc16_02, dims=1:16)
DimPlot(pbmc16_02, reduction = "tsne")
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

The method of choice for the visualization of the clusters is **“Uniform
manifold approximation and projection (UMAP)”.** It also tries to find a
low-dimensional representation that preserves relationships between
neighbors in high-dimensional space - that is the space made of the
number of principal components that have been used to find the
clusters - but also it produces more compact visual clusters with more
empty space between them. It also attempts to preserve more of the
global structure (distance) among clusters

``` r
pbmc16_02 <- RunUMAP(pbmc16_02, dims = 1:16)
DimPlot(pbmc16_02, reduction = "umap")
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

The presentation of the results is done in bidimensional space but the
computations were performed in the 16-dimensions, so using the
normalized and scaled expression values of each of the genes that make
the higher variability along the first 16 principal components.

What we can notice is that some clusters are very big, like the 0 and 1,
while the last two clusters are very small: this is something we expect
since from the PanglaoDB we know that there are some cell types for
which only few cells are present - like myocytes and Schwann cells. At
the same time these cell types are not reported as detected by the
article: the fact that they reported only 6 cell types may be due to the
fact that they decided to consider only those cell types that are
represented by a substantial number of cells coming from the muscle
sample.

**Number of cells for each cluster**

``` r
for (i in 0:11)
{print (paste('Cluster', i, 'has', length(Idents(pbmc16_02)[Idents(pbmc16_02) == i]), 'cells'))}
```

    ## [1] "Cluster 0 has 914 cells"
    ## [1] "Cluster 1 has 577 cells"
    ## [1] "Cluster 2 has 206 cells"
    ## [1] "Cluster 3 has 191 cells"
    ## [1] "Cluster 4 has 167 cells"
    ## [1] "Cluster 5 has 151 cells"
    ## [1] "Cluster 6 has 138 cells"
    ## [1] "Cluster 7 has 90 cells"
    ## [1] "Cluster 8 has 61 cells"
    ## [1] "Cluster 9 has 45 cells"
    ## [1] "Cluster 10 has 44 cells"
    ## [1] "Cluster 11 has 18 cells"

As expected, the number of cells in each cluster decreases moving from
cluster 0 to 11.

Now we can **check whether some of the critical quality parameters
influenced the clustering we got:** if we find a cluster with high
percentage of reads mapping on rbp genes or mitochondrial ones the risk
is that the cluster has as marker genes these types of genes that are
the less meaningful ones (it is a low quality cluster).

First we plot the distribution of the number of reads mapping in the
cells of each cluster

``` r
VlnPlot(pbmc16_02,features="nCount_RNA", pt.size=0)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

We can plot the distribution of the number of transcribed genes in the
cells of each cluster

``` r
VlnPlot(pbmc16_02,features="nFeature_RNA", pt.size=0)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

Plot the distribution of reads mapping on mitochondrial genes in the
cells of each cluster

``` r
VlnPlot(pbmc16_02,features="percent.mt", pt.size=0)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

We can see that cluster 9 has a very low percentage with respect to the
others, this can be sign of the fact that this cluster is created only
because its cells share the characteristic of having a particular low
number of reads mapping on mitochondrial genes so it an information to
take in account in the further assessment of the cell type of each
cluster.

Plot the distribution of reads mapping on ribosomal protein genes in the
cells of each cluster

``` r
VlnPlot(pbmc16_02,features="percent.rbp", pt.size=0)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

We can see that clusters 4 and 6 have a higher percentage of reads
mapping on rbp genes, this is something to take in consideration in the
further identification of their associated cell types in case we won’t
find significant marker genes.

Plot the distribution of cell cycle phase of the cells of each cluster

``` r
library(ggplot2)
pbmc16_02@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

The proportions of cells in the different cell cycle phases are quite
the same in all the clusters.

# Finding “marker” genes and assigning cell types to clusters

**Seurat includes a function that can be used to find genes
overexpressed in one clusters with respect to all the others.** To test
if a gene is overexpressed in one cluster there are different methods
but, for 10x data, the choice is to employ the **non parametric Wilcoxon
test.** We add another constraint to detect the marker genes: **to be a
marker gene the feature has to be expressed at least in 25% of the cells
of a cluster and with an avg_logFC threshold of at least 0.25.**

FindAllMarkers is an iterative function that, for each cluster, performs
the comparison with all the others

``` r
pbmc16_02.markers <- FindAllMarkers(pbmc16_02, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
```

Once the genes with the highest average expression in the cells of each
cluster (with respect to the cells of all the other clusters) are
identified, **we print for each cluster the first five more significant
markers sorted by avg_log2FC** - that is more significant than the
adjusted p-value because we are using a non parametric test that, when
run on very big clusters, may have the side effect of returning very low
adjusted p-values for small difference in numbers. So the advise is to
sort the putative marker genes by avg_log2FC

``` r
pbmc16_02.markers %>%group_by(cluster) %>%slice_max(n = 5, order_by = avg_log2FC)
```

    ## # A tibble: 60 × 7
    ## # Groups:   cluster [12]
    ##        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene   
    ##        <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>  
    ##  1 0               4.28 0.968 0.422 0         0       Fabp4  
    ##  2 5.73e-222       4.18 0.804 0.344 1.17e-217 0       Aqp1   
    ##  3 0               3.79 0.791 0.095 0         0       Gpihbp1
    ##  4 3.91e-238       3.65 0.757 0.198 7.96e-234 0       Id1    
    ##  5 0               3.26 0.832 0.113 0         0       Flt1   
    ##  6 1.77e-285       4.25 1     0.705 3.61e-281 1       Gsn    
    ##  7 0               4.00 0.981 0.101 0         1       Clec3b 
    ##  8 2.17e-303       3.91 0.827 0.097 4.42e-299 1       Smoc2  
    ##  9 0               3.76 0.99  0.22  0         1       Col3a1 
    ## 10 1.43e-157       3.72 0.53  0.074 2.92e-153 1       Ptx3   
    ## # … with 50 more rows

For each gene we have some quantities computed: **the p-value; the
adjusted p-value using BH adjustment; the avg_log2FC - that is the
average expression of each gene in the cells of a cluster divided by the
average expression in all the other cells not belonging to that
cluster - and the percentage of cells of the cluster expressing the gene
and the percentage of the cells outside the cluster expressing it.**

These are the genes that make the difference from one cluster to the
others and so it is important to search in which type of cell of the
muscle each gene is mostly expressed to associate to each cluster a cell
type. **All the genes have been found as markers of a specific cell type
on Tabula muris with exception of those indicated:**

| **Cluster** | **Marker**   | **Cell Type**                  |
|-------------|--------------|--------------------------------|
| 0           | Fabp4, Aqp1  | Endothelial cells              |
| 1           | Gsn, Clec3b  | Mesenchymal Stem cells         |
| 2           | Fmod, Thbs4  | Fibroblast \*                  |
| 3           | Sdc4, Myod1  | Skeletal muscle satellite cell |
| 4           | Cd79a, Cd74  | B cells                        |
| 5           | Lyz2, Retnla | Macrophage                     |
| 6           | Cd3g, Cd3d   | T cells                        |
| 7           | Rgs5, Acta2  | Smooth muscle cell \*\*        |
| 8           | Selp, Csf3   | Endothelial cells \*\*\*       |
| 9           | Acta1, Mylpf | Myocytes \*\*\*                |
| 10          | Nts, Lyve1   | Endothelial cells \*\*         |
| 11          | Mpz, Mbp     | Schwann cell \*\*\*            |

\*\* from The Human Protein Atlas: human homologs \*\*\* from PanglaoDB

**We can consider all together the top 10 marker genes - according to
the avg_log2FC - of each cluster and for each cell we plot its
expression value of all these genes in a heatmap**

``` r
pbmc16_02.markers %>%
group_by(cluster) %>%
top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc16_02, features = top10$gene, size=2) + theme(text = element_text(size = 3))+ NoLegend()
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->

In the heatmap the higher is the expression of a marker gene in a cell
the more the color moves to yellow. Looking at this plot we can see
which clusters tend to have in common some putative marker genes: for
example some cells of the cluster 1 seem to have a high expression also
of some marker genes of cluster 2 and 3. Anyway we can be satisfied of
the achieved clustering since **each cluster has a well defined yellow
rectangle meaning that a great percentage of cells of each cluster
express at the highest value only the marker genes of the cluster to
which they have been assigned.**

Since the marker genes of cluster 0 and 8 are those of endothelial cells
we can look for the markers common to cluster 0 and cluster 8 looking at
these two clusters against all the others

``` r
cluster0AND8.markers <- FindMarkers(pbmc16_02, ident.1 = c(0,8), min.pct = 0.25, test.use = "wilcox")
cluster0AND8.markers <- cluster0AND8.markers[order(-cluster0AND8.markers$avg_log2FC),]
head(cluster0AND8.markers, n = 10)
```

    ##                 p_val avg_log2FC pct.1 pct.2     p_val_adj
    ## Fabp4    0.000000e+00   5.430593 0.970 0.400  0.000000e+00
    ## Aqp1    1.327300e-248   4.985248 0.812 0.322 2.703976e-244
    ## Gpihbp1 1.621104e-307   3.841014 0.767 0.083 3.302513e-303
    ## Id1     4.138669e-230   3.665669 0.736 0.190 8.431297e-226
    ## Flt1     0.000000e+00   3.657963 0.841 0.080  0.000000e+00
    ## Cdh5     0.000000e+00   3.329694 0.811 0.082  0.000000e+00
    ## Esam     0.000000e+00   3.312679 0.823 0.104  0.000000e+00
    ## Cd36    8.317842e-301   3.240554 0.832 0.160 1.694511e-296
    ## Sox17   7.553853e-230   3.238760 0.600 0.045 1.538871e-225
    ## Rgcc    1.155788e-212   3.195707 0.799 0.341 2.354572e-208

**the top 10 marker genes, sorted by avg_logF2C, are the markers of the
endothelial cell type so we can conclude that cluster 0 and 8 are of the
same cell type: endothelial cells.**

Maybe these two clusters are subtypes of the same cell type: to
understand it we would look for genes overexpressed by cluster 8 but not
by cluster 0. Anyway endothelial cells don’t have specific subtypes so
we just merge the cells of these two clusters and it makes sense since
looking back at the UMAP we can see that cluster 8 is very close to
cluster 0.

In the same way, the marker genes of cluster 0 and 10 are those of
endothelial cells thus we can look for the markers common to cluster 0
and cluster 10 looking at these two clusters against all the others

``` r
cluster0AND10.markers <- FindMarkers(pbmc16_02, ident.1 = c(0,10), min.pct = 0.25, test.use = "wilcox")
cluster0AND10.markers <- cluster0AND10.markers[order(-cluster0AND10.markers$avg_log2FC),]
head(cluster0AND10.markers, n = 10)
```

    ##                 p_val avg_log2FC pct.1 pct.2     p_val_adj
    ## Fabp4    0.000000e+00   4.207386 0.938 0.425  0.000000e+00
    ## Aqp1    3.038305e-208   4.104381 0.785 0.343 6.189635e-204
    ## Gpihbp1 6.429118e-307   3.820328 0.770 0.088 1.309740e-302
    ## Id1     2.868716e-231   3.642586 0.741 0.193 5.844148e-227
    ## Cldn5   1.710787e-211   3.456370 0.594 0.066 3.485215e-207
    ## Flt1    4.919174e-295   3.181067 0.796 0.114 1.002134e-290
    ## Cdh5     0.000000e+00   3.133794 0.800 0.096  0.000000e+00
    ## Cav1    5.918091e-272   3.118586 0.867 0.340 1.205634e-267
    ## Timp4   1.740002e-213   3.107820 0.571 0.047 3.544732e-209
    ## Sox17   1.905291e-205   3.088064 0.582 0.061 3.881460e-201

**even in this case the top 10 marker genes, sorted by avg_logF2C, are
the markers of the endothelial cell type so we can conclude that cluster
0 and 10 are of the same cell type: endothelial cells.**

After these analysis we can merge cluster 8 and cluster 10 with cluster
0

``` r
Idents(pbmc16_02)[Idents(pbmc16_02) == 8] <- 0
Idents(pbmc16_02)[Idents(pbmc16_02) == 10] <- 0
head(Idents(pbmc16_02))
```

    ## AAACCTGAGCAGATCG AAACCTGGTCTAGTGT AAACCTGTCGCCTGAG AAACGGGCACACCGAC 
    ##                0                1                0                6 
    ## AAACGGGCATCTGGTA AAACGGGCATGTCTCC 
    ##                0                0 
    ## Levels: 0 1 2 3 4 5 6 7 9 11

We can now **plot the expression of some marker genes to see their
distribution in all the clusters**

``` r
VlnPlot(pbmc16_02, features = c("Fabp4", "Gsn"))
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

From the plot we can see that Fabp4 is expressed at high level only in
cells of cluster 0, while Gsn is expressed at the highest level only in
the cells of the cluster 1 and with a moderate value in some of the
cells of cluster 2.

``` r
DotPlot(pbmc16_02, features = c("Fabp4", "Lyz2","Sdc4", "Cd79a","Acta1","Gsn","Cd3g","Rgs5","Fmod","Mpz"))
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->

From the **dotplot** is visible that each of the plotted genes is a
unique marker of a cluster: in fact even if some genes, like Gsn, are
expressed in cells of other clusters too, each marker gene is expressed
with the highest average value and by the highest percentage of cells
only of one cluster.

We can **plot a UMAP for each marker gene** to see which cluster of
cells express each marker

``` r
FeaturePlot(pbmc16_02, features = c("Fabp4", "Gsn", "Fmod","Sdc4", "Cd79a", "Lyz2","Cd3g","Rgs5","Acta1","Mpz"))
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-54-1.png)<!-- -->

Almost all the markers are expressed only in one cluster of cells with
the exception of Gsn and Sdc4 that - as already seen in the dotplot -
even if expressed with highest value by all the cells of one cluster are
also highly expressed by a smaller percentage of cells of other
clusters.

Eventually we can **plot the cells with the corresponding cell type
label**

``` r
new.cluster.ids <- c("endothelial cell", "mesenchymal stem cell", "fibroblast", "skeletal muscle satellite cell", "B", "macrophage", "T", "smooth muscle cell", "myocytes", "Schwann cell")
names(new.cluster.ids) <- levels(pbmc16_02)
pbmc16_02 <- RenameIdents(pbmc16_02, new.cluster.ids)
DimPlot(pbmc16_02, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-55-1.png)<!-- -->

We **were able to identify all the cell types suggested by PanglaoDB but
two - Pericytes and Myoblasts. We detected three cell types -
Mesenchymal stem cells, Skeletal muscle satellite cell and Smooth muscle
cell - that they don’t report** (while they report some unknown cells
that we were able to associate with specific cell types). This
difference is probably due to the fact that they performed the analysis
with some “standard” parameters that may be not the best one for this
data.

We **identified all the cell types detected by the article plus other
four cell types: Fibroblasts, Schwann cells, Smooth muscle cell and
Myocytes.** This difference with the article analysis is due to the fact
that they started from a different dataset than ours.

# Same analysis with different parameters

Starting from the same pbmc Seurat object - previously filtered to
remove all the cells with a number of features smaller than 200 or
higher than 4000 as well as those with more than 4% of reads mapping on
mitochondrial genes - we now **set new parameters for the clustering to
see if and how the final results changes:**

-   Neighboring cells are computed along the first 25 principal
    components, which, as seen from the previous Elbow plot, corresponds
    to the point from which the curve starts to flatten and the SD is
    lower than 2;
-   Clustering Resolution is set to 0.01 to reduce the number of
    clusters.

``` r
pbmc25_001 <- FindNeighbors(pbmc, dims = 1:25)
pbmc25_001 <- FindClusters(pbmc25_001, resolution = 0.01)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2602
    ## Number of edges: 83893
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9954
    ## Number of communities: 7
    ## Elapsed time: 0 seconds

We can look at the cluster IDs of the first 5 cells

``` r
head(Idents(pbmc25_001), 5)
```

    ## AAACCTGAGCAGATCG AAACCTGGTCTAGTGT AAACCTGTCGCCTGAG AAACGGGCACACCGAC 
    ##                0                1                0                2 
    ## AAACGGGCATCTGGTA 
    ##                0 
    ## Levels: 0 1 2 3 4 5 6

With these parameters **only 7 clusters are identified.**

**Clusters can be visualized using tSNE**

``` r
pbmc25_001 <- RunTSNE(pbmc25_001, dims=1:25)
DimPlot(pbmc25_001, reduction = "tsne")
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->

**Clusters can be visualized using UMAP**

``` r
pbmc25_001 <- RunUMAP(pbmc25_001, dims = 1:25)
DimPlot(pbmc25_001, reduction = "umap")
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-59-1.png)<!-- -->

**Number of cells for each cluster**

``` r
for (i in 0:6)
{print (paste('Cluster', i, 'has', length(Idents(pbmc25_001)[Idents(pbmc25_001) == i]), 'cells'))}
```

    ## [1] "Cluster 0 has 1000 cells"
    ## [1] "Cluster 1 has 796 cells"
    ## [1] "Cluster 2 has 453 cells"
    ## [1] "Cluster 3 has 191 cells"
    ## [1] "Cluster 4 has 89 cells"
    ## [1] "Cluster 5 has 45 cells"
    ## [1] "Cluster 6 has 28 cells"

# Overlapping between the two results

We recompute the clusters as they are computed looking only at the first
16 PCs and with a resolution parameters of 0.2

``` r
pbmc16_02 <- FindNeighbors(pbmc, dims = 1:16)
pbmc16_02 <- FindClusters(pbmc16_02, resolution = 0.2)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2602
    ## Number of edges: 80813
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9586
    ## Number of communities: 12
    ## Elapsed time: 0 seconds

Then **we compare the two results obtained using different parameters:
we compute how many cells - each identified by its barcode - are shared
by each possible pair of clusters**

``` r
overlap <- data.frame(cluster_res2 = NA, cluster_res1 = NA, number_cells_common =NA)

for (i in 0:6) {
  for (j in 0:11) {
    overlap <- rbind(overlap, c(i,j,length(intersect(rownames(as.data.frame(Idents(pbmc16_02)[Idents(pbmc16_02) == j])),rownames(as.data.frame(Idents(pbmc25_001)[Idents(pbmc25_001) == i]))))))
  }
}

overlap <- overlap[c(-1),]
head(overlap)
```

    ##   cluster_res2 cluster_res1 number_cells_common
    ## 2            0            0                 914
    ## 3            0            1                   0
    ## 4            0            2                   2
    ## 5            0            3                   0
    ## 6            0            4                   0
    ## 7            0            5                   3

We can **visualize through a dotplot the result of the comparison**

``` r
ggplot(overlap, aes(x= cluster_res2, y= cluster_res1, colour=number_cells_common, size=number_cells_common)) +  geom_point() +
  theme_classic() +  scale_x_continuous(breaks = 0:6) + scale_y_continuous(breaks = 0:11)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-63-1.png)<!-- -->

Cluster 0 of result 2 contains all the cells of clusters 0, 8 and 11 of
result 1. Cluster 1 of result 2 contains all the cells of clusters 1 and
2 of result 1 and part of the cells of cluster 10. Cluster 2 of result 2
contains all the cells of clusters 4, 5 and 6 of result 1. Cluster 3 of
result 2 contains all the cells of clusters 3 of result 1. Cluster 4 of
result 2 contains all the cells of cluster 7 of result 1. Cluster 5 of
result 2 contains all the cells of cluster 9 of result 1. Cluster 6 of
result 2 contains part of the cells of cluster 10 of result 1.

We should now consider the same QC parameters as before to observe
whether any of the clusters is of low quality.

Distribution of reads per cell per each cluster:

``` r
VlnPlot(pbmc25_001,features="nCount_RNA", pt.size=0)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-64-1.png)<!-- -->

Distribution of transcribed genes per each cluster:

``` r
VlnPlot(pbmc25_001,features="nFeature_RNA", pt.size=0)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-65-1.png)<!-- -->

Distribution of the percentage of reads mapping on mitochondrial genes
per cluster:

``` r
VlnPlot(pbmc25_001,features="percent.mt", pt.size=0)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-66-1.png)<!-- -->

We can see that cluster 5 has a very low percentage with respect to the
others, this can be sign of the fact that this cluster is created only
because its cells share the characteristic of having a particular low
number of reads mapping on mitochondrial genes so it an information to
take in account in the further assessment of the cell type of each
cluster.

Distribution of the percentage of reads mapping on rRNA genes per
cluster:

``` r
VlnPlot(pbmc25_001,features="percent.rbp",pt.size=0)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-67-1.png)<!-- -->

Notice the higher percentage of ribosomal genes in clusters 2 and 3
compared to the others. In case the marker genes for these clusters are
indeed ribosomal protein genes, the clusters may be deemed of low
quality and will be discarded.

Display the portions of cells per cell cycle phase in each cluster:

``` r
library(ggplot2)
pbmc25_001@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-68-1.png)<!-- -->

Cells in cluster 5 appear not to be proliferating as much as those in
the other clusters.

We can now **retrieve the marker genes for each cluster.** These are
selected as the 5 genes having the highest average Log_2 Fold Change of
expression in each cluster

``` r
pbmc25_001.markers <- FindAllMarkers(pbmc25_001, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc25_001.markers %>%group_by(cluster) %>%slice_max(n = 5, order_by = avg_log2FC)
```

    ## # A tibble: 35 × 7
    ## # Groups:   cluster [7]
    ##        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene   
    ##        <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>  
    ##  1 0               5.42 0.952 0.403 0         0       Fabp4  
    ##  2 4.27e-235       4.94 0.795 0.325 8.69e-231 0       Aqp1   
    ##  3 2.94e-293       3.80 0.748 0.084 5.99e-289 0       Gpihbp1
    ##  4 7.70e-220       3.63 0.721 0.191 1.57e-215 0       Id1    
    ##  5 0               3.62 0.821 0.081 0         0       Flt1   
    ##  6 0               5.45 0.997 0.395 0         1       Dcn    
    ##  7 0               5.09 1     0.669 0         1       Gsn    
    ##  8 0               4.65 0.979 0.121 0         1       Col1a1 
    ##  9 0               4.46 0.984 0.159 0         1       Col1a2 
    ## 10 3.11e- 61       4.42 0.474 0.182 6.33e- 57 1       Fmod   
    ## # … with 25 more rows

Annotation of cell types is done manually by looking up the marker genes
in the Tabula Muris database. In the following table, only two out of
the five marker genes are reported.

| **Cluster** | **Marker**   | **Cell Type**                  |
|-------------|--------------|--------------------------------|
| 0           | Fabp4, Aqp1  | Endothelial cells              |
| 1           | Dcn,Gsn      | Mesenchymal Stem cells         |
| 2           | Cd74, H2-Aa  | B cells                        |
| 3           | Sdc4, Myod1  | Skeletal muscle satellite cell |
| 4           | Rgs5, Acta2  | Vascular smooth muscle cells   |
| 5           | Acta1, Mylpf | Skeletal myocytes \*\*         |
| 6           | Nts, Lyve1   | Endothelial cells \*\*\*       |

\*\* identification by querying Google search engine ‘Gene is expressed
in xxx mouse cells’ \*\*\* from The Human Protein Atlas: human homologs

Show a **heatmap of the expression of the top 10 marker genes for each
cluster**

``` r
pbmc25_001.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc25_001, features = top10$gene, size=2) + theme(text = element_text(size = 3))+ NoLegend()
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-70-1.png)<!-- -->

The overall clustering seems satisfactory as the well defined yellow
rectangles indicate that genes are indeed specific for each cluster.

We can also use a **dotplot to observe the expression of some marker
genes with respect to clusters:**

``` r
DotPlot(pbmc25_001, features = c("Fabp4", "Dcn","Cd74", "Sdc4","Rgs5","Acta1","Nts"))
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-71-1.png)<!-- -->

Nevertheless, cluster 6 and 0 both seem to host endothelial cells based
on the manual annotation performed, we therefore check whether it is the
case by firstly checking which are the markers of the two clusters
considered jointly:

``` r
cluster0AND6.markers <- FindMarkers(pbmc25_001, ident.1 = c(0,6), min.pct = 0.25, test.use = "wilcox")
cluster0AND6.markers <- cluster0AND6.markers[order(-cluster0AND6.markers$avg_log2FC),]
head(cluster0AND6.markers, n = 10)
```

    ##                 p_val avg_log2FC pct.1 pct.2     p_val_adj
    ## Fabp4    0.000000e+00   5.396023 0.937 0.403  0.000000e+00
    ## Aqp1    6.083720e-228   4.911712 0.785 0.323 1.239375e-223
    ## Gpihbp1 3.174857e-295   3.878685 0.742 0.076 6.467818e-291
    ## Id1     1.221109e-219   3.639948 0.716 0.185 2.487642e-215
    ## Flt1     0.000000e+00   3.582841 0.801 0.081  0.000000e+00
    ## Cldn5   1.748556e-190   3.428840 0.558 0.065 3.562158e-186
    ## Cdh5     0.000000e+00   3.381375 0.793 0.069  0.000000e+00
    ## Esam    5.600750e-307   3.260481 0.793 0.099 1.140985e-302
    ## Pecam1   0.000000e+00   3.235256 0.779 0.081  0.000000e+00
    ## Egfl7   3.048999e-272   3.233692 0.696 0.064 6.211421e-268

As the best marker genes are specific of endothelial cell, we can be
sure that both clusters are of the same cell type. We now merge the two
clusters:

``` r
Idents(pbmc25_001)[Idents(pbmc25_001) == 6] <- 0
head(Idents(pbmc25_001))
```

    ## AAACCTGAGCAGATCG AAACCTGGTCTAGTGT AAACCTGTCGCCTGAG AAACGGGCACACCGAC 
    ##                0                1                0                2 
    ## AAACGGGCATCTGGTA AAACGGGCATGTCTCC 
    ##                0                0 
    ## Levels: 0 1 2 3 4 5

Plot the resulting cell types over the UMAP representation:

``` r
new.cluster.ids <- c("endothelial cell", "mesenchymal stem cell", "B", "skeletal muscle satellite cell", "vascular smooth muscle cell", "skeletal myocytes")
names(new.cluster.ids) <- levels(pbmc25_001)
pbmc25_001 <- RenameIdents(pbmc25_001, new.cluster.ids)
DimPlot(pbmc25_001, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-74-1.png)<!-- -->

We can see that the **result of the clustering performed considering the
first 25 PCs and a resolution parameter of 0.01 is for sure worse than
the one obtained using PC=16 and resolution 0.2 since:**

-   from the UMAP is visible that the cells of some cell types - as B
    cells and mesenchymal stem cells - are divided in more than one
    cluster;
-   some of the cell types that have been previously identified - like T
    cells, Schwann cells, macrophages and fibroblasts - are not
    highlighted by this solution meaning that the clustering performed
    with these parameters is not able to distinguish among some of the
    cell types.

In fact from the overlapping of the clusters obtained with the two
analysis with different parameters it results that:

-   all Schwann cells - that are those of cluster 11 in the first
    analysis - are considered endothelial cells by the second one;
-   all Fibroblasts - that are those of cluster 2 in the first
    analysis - are considered Mesenchymal Stem cells by the second one
    (that is the reason why in the UMAP there are two clusters
    associated to Mesenchymal stem cells);
-   all Macrophages and T cells - that are those of cluster 5 and 6 in
    the first analysis - are considered B cells by the second one (that
    is the reason why in the UMAP there are three clusters associated to
    B cells).

**In conclusion these parameters, PC=25 and resolution 0.01, are not
good since they allow to identify a reduced number of clusters in which
the cells are not all of the same cell type.**

# Same analysis performed using the article parameters

In the article from which this data comes from some modifications to the
previous pipeline and parameters were introduced, so we decided to
repeat the analysis applying these modifications.

First we clean the workspace and load the count table

``` r
rm(list=ls())
load('SRA653146_SRS3044259.sparse.RData')
```

Cut the name of each gene as before

``` r
rownames(sm) <- lapply(rownames(sm), FUN = function(x) {
      if (str_count(x, "_") != 1) {
        a <- strsplit(x, '_')
        x <- paste(a[[1]][1],a[[1]][2],sep='_')
  }
  else {
        x <-strsplit(x, '_')[[1]][1]}
  })
  
head(rownames(sm))
```

    ## [1] "00R_AC107638.2" "0610009O20Rik"  "1010001N08Rik"  "1110020A21Rik" 
    ## [5] "1700007L15Rik"  "1700015F17Rik"

Initialize the Seurat object with the correct gene symbols

``` r
pbmc <- CreateSeuratObject(
     sm,
     project = "scRNA_seq_dermis",
     min.cells = 3,
     min.features = 200
)
```

Add inside the Seurat object the percentage of reads mapping on
mitochondrial and rbp genes

``` r
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
pbmc[["percent.rbp"]] <- PercentageFeatureSet(pbmc, pattern = "^Rp[Ls]")
```

## Quality control

We plot the distribution of the different quantities

``` r
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4, pt.size=0)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-79-1.png)<!-- -->

``` r
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-80-1.png)<!-- -->

``` r
plot2
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-81-1.png)<!-- -->

As seen before there is no correlation between the percentage of reads
mapping on mitochondrial genes and the total number of reads from each
cell. Instead there is a strong positive correlation between the number
of reads coming from each cells and the amount of its transcribed genes.

``` r
plot3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.rbp")
plot3
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-82-1.png)<!-- -->

As before the is no correlation between the percentage of reads mapping
on rbp genes and the total number of reads from each cell.

**As before we use the same thresholds to filter our data, furthermore
the article suggests to remove the cells with less than 500 expressed
genes and less than 500 reads**

``` r
pbmc <- subset(pbmc, subset = ((nFeature_RNA > 500 & nFeature_RNA < 4000 ) & nCount_RNA > 500 & percent.mt < 4))
```

Let us see how many samples remain:

``` r
pbmc
```

    ## An object of class Seurat 
    ## 20372 features across 2447 samples within 1 assay 
    ## Active assay: RNA (20372 features, 0 variable features)

Now more cells are thrown away.

## Normalizing the data

From the github code that is indicated in the article
(<https://github.com/czbiohub/tabula-muris/blob/master/00_data_ingest/02_tissue_analysis_rmd/Limb_Muscle_facs.Rmd>)
we can see that **first they normalized the counts, then they scaled and
eventually selected the most variable genes but in a different way than
what was previously done**

``` r
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```

Associate to each cell its cell cycle phase

``` r
CellCycleScoring(pbmc, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> pbmc
```

**Scaling of the normalized counts**

``` r
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

**The selection of the most variable genes is performed in this way by
the authors**

``` r
pbmc <- FindVariableFeatures(pbmc, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5)
```

**Run PCA** to find the most variable genes in each principal component

``` r
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```

    ## PC_ 1 
    ## Positive:  Col6a1, Dcn, Serping1, Col1a2, Gsn, Col1a1, Igfbp6, Col3a1, Col6a2, Pcolce 
    ##     Serpinf1, Mmp2, Fstl1, Timp2, Clec3b, Rnase4, Mfap5, Ogn, Nbl1, Il11ra1 
    ##     Col5a2, Dpt, Axl, Islr, Htra3, Bgn, Lum, Col6a3, Fndc1, Tnxb 
    ## Negative:  Srgn, Tmsb4x, Cd52, Fabp4, Coro1a, Laptm5, Ptpn18, Rac2, Cd53, Gpihbp1 
    ##     Gmfg, Aqp1, H2-Aa, Psmb8, Cd37, Cd74, Ctss, H2-Eb1, H2-Ab1, Arhgap45 
    ##     Tinagl1, Arhgdib, Egfl7, H2-DMa, Emcn, Cytip, C1qtnf9, RP23-310J6.1, Ppp1r16b, Ctla2a 
    ## PC_ 2 
    ## Positive:  Fabp4, Hspb1, Tinagl1, Gpihbp1, Id1, Aqp1, Apold1, Egfl7, Emcn, RP23-310J6.1 
    ##     Tm4sf1, Rgcc, C1qtnf9, Sox17, Rasip1, Timp4, Cxcl12, Tcf15, Rnd1, Cldn5 
    ##     Iigp1, Ctla2a, Fabp5, Grrp1, Icam2, Isg15, Id3, Cdc42ep3, Tsc22d1, Cd93 
    ## Negative:  Ctss, Tyrobp, Fcer1g, Alox5ap, Lyz2, Laptm5, Cd74, H2-Ab1, H2-Eb1, H2-DMa 
    ##     H2-Aa, Ccl9, Spi1, Fcgr3, Cd53, Cd52, Cd68, H2-DMb1, C1qb, C1qc 
    ##     Ccl6, Ly86, Ptafr, Csf1r, Wfdc17, Cd83, Pld4, C1qa, Ctsc, Lilr4b 
    ## PC_ 3 
    ## Positive:  Cst3, C1qa, C1qb, C1qc, Csf1r, Lyz2, Fcgr3, Fcer1g, Pf4, Ccl9 
    ##     Alox5ap, Mrc1, C5ar1, Cd68, Ccl6, Cd14, Ptafr, C3ar1, Ier3, F13a1 
    ##     Fcrls, Tyrobp, Fabp4, Lilr4b, Adgre1, Cbr2, Aif1, Cfp, Mgl2, P2ry6 
    ## Negative:  Satb1, Ptprcap, RP24-146B4.2, Cd79a, Ltb, Ccr7, Limd2, RP23-32C18.6, Cd79b, Ighm 
    ##     Rac2, Ighd, Igkc, Arhgap45, Ly6d, Fcmr, RP24-547D11.2, Iglc2, Cd2, H2-DMb2 
    ##     Coro1a, Gimap3, Trbc2, Cd3d, Cd52, H2-Ob, Cd3g, Ms4a1, RP24-338A5.4, H2-Oa 
    ## PC_ 4 
    ## Positive:  Cilp2, Tnmd, Fmod, Chad, Thbs4, Comp, Col11a1, Cpxm2, Ptx4, Col12a1 
    ##     Kera, Itgbl1, Mfap4, Col8a2, AC159092.1, Scx, Col11a2, Angptl7, Tnc, Wif1 
    ##     Cilp, Mkx, Ptgis, Lox, Olfml3, Clec11a, Wisp1, Kctd1, Fxyd6, Timp1 
    ## Negative:  Dpt, Cd34, Ifi205, Clec3b, Scara5, Cd248, Ccl11, Dpep1, Efemp1, Cadm3 
    ##     Ifi204, Mfap5, Plac9c, Pi16, Emilin2, Efhd1, Mnda, Ace, Pla1a, Ptx3 
    ##     Entpd2, Col14a1, Pdgfra, Ecm1, Ifi207, Cd55, Pcsk6, Gas7, Nid1, Ms4a4d 
    ## PC_ 5 
    ## Positive:  Tmsb4x, Actb, Tnmd, Ptprcap, Cilp2, Fmod, Tubb5, Limd2, Ltb, Satb1 
    ##     Thbs4, RP24-146B4.2, Itgbl1, Arhgap45, Ccr7, Rac2, Cpxm2, Coro1a, Srgn, Chad 
    ##     Cilp, Comp, Emb, Col12a1, Cd52, Mfap4, Kera, Ptx4, Sh3bgrl3, Ptgis 
    ## Negative:  Myl1, Ckm, RP24-426K19.3, Pgam2, Cox6a2, Atp2a1, Tcap, Smpx, Tnni2, Cox8b 
    ##     Eno3, Acta1, Myoz1, Tnnt3, Pvalb, Mylpf, Eef1a2, Car3, Actn3, Apobec2 
    ##     Adssl1, Tmod4, Rpl3l, Myot, Sh3bgr, Cox7a1, Mybpc2, Ak1, Tpm2, Ppp1r27

``` r
print(pbmc[["pca"]], dims = 1:5, , nfeatures = 5)
```

    ## PC_ 1 
    ## Positive:  Col6a1, Dcn, Serping1, Col1a2, Gsn 
    ## Negative:  Srgn, Tmsb4x, Cd52, Fabp4, Coro1a 
    ## PC_ 2 
    ## Positive:  Fabp4, Hspb1, Tinagl1, Gpihbp1, Id1 
    ## Negative:  Ctss, Tyrobp, Fcer1g, Alox5ap, Lyz2 
    ## PC_ 3 
    ## Positive:  Cst3, C1qa, C1qb, C1qc, Csf1r 
    ## Negative:  Satb1, Ptprcap, RP24-146B4.2, Cd79a, Ltb 
    ## PC_ 4 
    ## Positive:  Cilp2, Tnmd, Fmod, Chad, Thbs4 
    ## Negative:  Dpt, Cd34, Ifi205, Clec3b, Scara5 
    ## PC_ 5 
    ## Positive:  Tmsb4x, Actb, Tnmd, Ptprcap, Cilp2 
    ## Negative:  Myl1, Ckm, RP24-426K19.3, Pgam2, Cox6a2

**Plot the projection of the cells in the first two principal
components**

``` r
DimPlot(pbmc, reduction = "pca")
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-91-1.png)<!-- -->

We can see that cells in different CC phases are all mixed so there
isn’t a cell cycle effect.

## Clustering

``` r
ElbowPlot(pbmc, ndims=50)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-92-1.png)<!-- -->

As before, from the Elbow plot we can say that the best number of PC
parameter is 25/30 since for this value the SD of the data that is not
explained is less than 20% and it is the starting point of the plateau.
But in **the article it is suggested to consider the first 10 PCs and so
we will now test this parameter.**

We **find the k neighbor cells** of each cell

``` r
pbmc10_05 <- FindNeighbors(pbmc, dims = 1:10)
```

**The Louvain algorithm is used to perform clustering with a resolution
of 0.5 (suggested by the article)**

``` r
pbmc10_05 <- FindClusters(pbmc10_05, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2447
    ## Number of edges: 76972
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9211
    ## Number of communities: 12
    ## Elapsed time: 0 seconds

``` r
head(Idents(pbmc10_05), 5)
```

    ## AAACCTGAGCAGATCG AAACCTGGTCTAGTGT AAACCTGTCGCCTGAG AAACGGGCACACCGAC 
    ##                6                1                0                9 
    ## AAACGGGCATCTGGTA 
    ##                0 
    ## Levels: 0 1 2 3 4 5 6 7 8 9 10 11

**12 clusters are identified as in the analysis done with PC=16 and
resolution of 0.2.**

``` r
DimPlot(pbmc10_05, reduction = "pca")
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-96-1.png)<!-- -->

``` r
pbmc10_05 <- RunTSNE(pbmc10_05, dims=1:10)
DimPlot(pbmc10_05, reduction = "tsne")
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-97-1.png)<!-- -->

``` r
pbmc10_05 <- RunUMAP(pbmc10_05, dims = 1:10)
DimPlot(pbmc10_05, reduction = "umap")
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-98-1.png)<!-- -->

**The clusters are not so well separated one from the others as before
and different clusters, like 1 - 5 or 0 - 2 - 6 are very close.**

``` r
for (i in 0:11)
{print (paste('Cluster', i, 'has', length(Idents(pbmc10_05)[Idents(pbmc10_05) == i]), 'cells'))}
```

    ## [1] "Cluster 0 has 476 cells"
    ## [1] "Cluster 1 has 393 cells"
    ## [1] "Cluster 2 has 281 cells"
    ## [1] "Cluster 3 has 200 cells"
    ## [1] "Cluster 4 has 191 cells"
    ## [1] "Cluster 5 has 175 cells"
    ## [1] "Cluster 6 has 172 cells"
    ## [1] "Cluster 7 has 167 cells"
    ## [1] "Cluster 8 has 145 cells"
    ## [1] "Cluster 9 has 137 cells"
    ## [1] "Cluster 10 has 89 cells"
    ## [1] "Cluster 11 has 21 cells"

``` r
VlnPlot(pbmc10_05,features="nCount_RNA", pt.size = 0)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-100-1.png)<!-- -->

``` r
VlnPlot(pbmc10_05,features="nFeature_RNA", pt.size = 0)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-101-1.png)<!-- -->

``` r
VlnPlot(pbmc10_05,features="percent.mt", pt.size = 0)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-102-1.png)<!-- -->

The percentage of reads mapping on mitochondrial genes is low in all the
clusters.

``` r
VlnPlot(pbmc10_05,features="percent.rbp", pt.size = 0)
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-103-1.png)<!-- -->

We can see that cluster 4,7 and 9 have a higher percentage of reads
mapping on rbp genes, this is something to take in consideration in the
further identification of their associated cell type in case we won’t
find significant marker genes.

``` r
library(ggplot2)
pbmc10_05@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-104-1.png)<!-- -->

The proportion of cells in different CC phases is more or less the same
in all the clusters.

## Finding “marker” genes and assigning cell types to clusters

``` r
pbmc10_05.markers <- FindAllMarkers(pbmc10_05, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
```

``` r
pbmc10_05.markers %>%group_by(cluster) %>%slice_max(n = 5, order_by = avg_log2FC)
```

    ## # A tibble: 60 × 7
    ## # Groups:   cluster [12]
    ##        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene        
    ##        <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>       
    ##  1 6.82e-208       2.79 0.79  0.155 1.39e-203 0       Tcf15       
    ##  2 3.68e-205       2.73 0.832 0.195 7.51e-201 0       Egfl7       
    ##  3 6.46e-159       2.68 0.889 0.459 1.32e-154 0       Sgk1        
    ##  4 9.29e-173       2.66 0.828 0.246 1.89e-168 0       RP23-310J6.1
    ##  5 3.41e-175       2.65 0.987 0.778 6.94e-171 0       Klf2        
    ##  6 8.15e-292       4.13 0.949 0.144 1.66e-287 1       Smoc2       
    ##  7 1.11e-122       3.66 0.621 0.134 2.26e-118 1       Cxcl14      
    ##  8 3.80e-194       3.31 1     0.762 7.73e-190 1       Gsn         
    ##  9 2.39e-239       3.15 0.916 0.168 4.87e-235 1       Lum         
    ## 10 1.20e-196       2.72 0.858 0.163 2.44e-192 1       Myoc        
    ## # … with 50 more rows

Looking at the top 5 expressed genes in each cluster we can see that for
some clusters they are different from the previous ones and in fact the
result is quite different, in particular:

| **Cluster** | **Marker**    | **Cell Type**                  |
|-------------|---------------|--------------------------------|
| 0           | Tcf15, Egfl7  | Endothelial cells              |
| 1           | Smoc2, Cxcl14 | Mesenchymal Stem cells         |
| 2           | Ifit1, Iigp1  | Endothelial cells              |
| 3           | Fmod, Thbs4   | Fibroblasts                    |
| 4           | Sdc4, Myod1   | Skeletal muscle satellite cell |
| 5           | Pi16, Mfap5v  | Mesenchymal stem cells         |
| 6           | Selp, Csf3    | Endothelial cells              |
| 7           | Cd79a, Cd74   | B cells                        |
| 8           | Lyz2, Retnla  | Macrophage                     |
| 9           | Cd3g, Cd3d    | T cells                        |
| 10          | Rgs5, Acta2   | Smooth muscle cell \*\*        |
| 11          | Acta1, Tnnt3  | Myocytes \*\*\*                |

\*\* from The Human Protein Atlas: human homologs \*\*\* from PanglaoDB

``` r
pbmc10_05.markers %>%
group_by(cluster) %>%
top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc10_05, features = top10$gene, size=2) + theme(text = element_text(size = 3))+ NoLegend()
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-107-1.png)<!-- -->

The heatmap shows more cases than before in which some cells of a
cluster have a high expression also of some marker genes of other
clusters: for example some cells of cluster 0 highly express marker
genes of cluster 2, while some cells of cluster 2 express marker genes
of cluster 0 and 6 and some of cluster 5 express marker genes of
cluster 1. The same between cluster 2 and cluster 6. Anyway **the
overall result is quite good but not as good as the one obtained using
PC=16 and resolution of 0.2.**

Since the marker genes of cluster 0, 2 and 6 are those of endothelial
cells we can find the markers common to cluster 0 and cluster 2 looking
at cluster 0 and 2 against all the other clusters

``` r
cluster0AND2.markers <- FindMarkers(pbmc10_05, ident.1 = c(0,2), min.pct = 0.25, test.use = "wilcox")
cluster0AND2.markers <- cluster0AND2.markers[order(-cluster0AND2.markers$avg_log2FC),]
head(cluster0AND2.markers, n = 10)
```

    ##                      p_val avg_log2FC pct.1 pct.2     p_val_adj
    ## Aqp1         3.174770e-238   4.271760 0.859 0.352 6.467642e-234
    ## Fabp4         0.000000e+00   4.121901 0.999 0.440  0.000000e+00
    ## Gpihbp1       0.000000e+00   3.859615 0.882 0.098  0.000000e+00
    ## Id1          6.989438e-255   3.617443 0.827 0.213 1.423888e-250
    ## Rgcc         3.048920e-259   3.246073 0.910 0.352 6.211260e-255
    ## Timp4        4.947658e-261   3.245739 0.678 0.051 1.007937e-256
    ## Cd36          0.000000e+00   3.181412 0.914 0.196  0.000000e+00
    ## Flt1          0.000000e+00   3.156052 0.906 0.131  0.000000e+00
    ## RP23-310J6.1 2.315692e-234   3.128455 0.778 0.172 4.717528e-230
    ## Cav1         9.147214e-284   3.054203 0.937 0.374 1.863470e-279

**the top 10 marker genes, sorted by avg_logF2C, are the markers of the
endothelial cell type so cluster 0 and 2 are of the same cell type:
endothelial cells.** In fact cluster 2 is close to cluster 0 looking at
UMAP.

In the same way, we can find the markers common to cluster 0 and cluster
6 looking at cluster 0 and 6 against all the other clusters

``` r
cluster0AND6.markers <- FindMarkers(pbmc10_05, ident.1 = c(0,6), min.pct = 0.25, test.use = "wilcox")
cluster0AND6.markers <- cluster0AND6.markers[order(-cluster0AND6.markers$avg_log2FC),]
head(cluster0AND6.markers, n = 10)
```

    ##                      p_val avg_log2FC pct.1 pct.2     p_val_adj
    ## Egfl7        6.562284e-219   2.700566 0.775 0.155 1.336868e-214
    ## Tcf15        1.812072e-148   2.411866 0.637 0.150 3.691554e-144
    ## Fabp4        1.189107e-172   2.391010 0.915 0.504 2.422448e-168
    ## Sgk1         3.420081e-113   2.332146 0.778 0.457 6.967389e-109
    ## Klf2         6.739341e-128   2.308664 0.935 0.777 1.372939e-123
    ## Cav1         1.394156e-180   2.300935 0.900 0.421 2.840174e-176
    ## Cav2         2.735454e-181   2.286986 0.789 0.228 5.572666e-177
    ## Sdpr         9.083495e-195   2.266614 0.861 0.263 1.850490e-190
    ## Gpihbp1      9.138376e-162   2.264635 0.744 0.195 1.861670e-157
    ## RP23-310J6.1 1.091120e-125   2.260132 0.696 0.238 2.222830e-121

**the top 10 marker genes, sorted by avg_logF2C, are the markers of the
endothelial cell type so cluster 0 and 6 are of the same cell type:
endothelial cells.** In fact even cluster 6 is close to cluster 0
looking at UMAP.

Since the marker genes of cluster 1 and 5 are those of Mesenchymal stem
cell we can find the markers common to cluster 1 and cluster 5 looking
at cluster 1 and 5 against all the other clusters

``` r
cluster1AND5.markers <- FindMarkers(pbmc10_05, ident.1 = c(1,5), min.pct = 0.25, test.use = "wilcox")
cluster1AND5.markers <- cluster1AND5.markers[order(-cluster1AND5.markers$avg_log2FC),]
head(cluster1AND5.markers, n = 10)
```

    ##                p_val avg_log2FC pct.1 pct.2     p_val_adj
    ## Gsn    3.914904e-276   4.219734 1.000 0.740 7.975443e-272
    ## Clec3b  0.000000e+00   3.985462 0.984 0.109  0.000000e+00
    ## Smoc2  1.109810e-284   3.900056 0.831 0.104 2.260904e-280
    ## Ptx3   1.168422e-149   3.777839 0.535 0.078 2.380309e-145
    ## Col3a1  0.000000e+00   3.740890 0.989 0.235  0.000000e+00
    ## Lum     0.000000e+00   3.360450 0.901 0.103  0.000000e+00
    ## Ifi205  0.000000e+00   3.343752 0.820 0.070  0.000000e+00
    ## Pi16   1.958198e-160   3.312820 0.613 0.104 3.989241e-156
    ## Dpt     0.000000e+00   3.294717 0.938 0.117  0.000000e+00
    ## Myoc   1.868717e-272   3.166370 0.826 0.108 3.806951e-268

\*\*the top 10 marker genes, sorted by avg_logF2C, are the markers of
Mesenchymal stem cell type so cluster 1 and 5 are of the same cell type:
endothelial cells.\* In fact cluster 5 is close to cluster 1 looking at
UMAP.

After these analysis we can merge the cells of cluster 0, 2 and 6 and
also cells of cluster 1 and 5

``` r
Idents(pbmc10_05)[Idents(pbmc10_05) == 2] <- 0
Idents(pbmc10_05)[Idents(pbmc10_05) == 6] <- 0
Idents(pbmc10_05)[Idents(pbmc10_05) == 5] <- 1
head(Idents(pbmc10_05))
```

    ## AAACCTGAGCAGATCG AAACCTGGTCTAGTGT AAACCTGTCGCCTGAG AAACGGGCACACCGAC 
    ##                0                1                0                9 
    ## AAACGGGCATCTGGTA AAACGGGGTGTCTGAT 
    ##                0                0 
    ## Levels: 0 1 3 4 7 8 9 10 11

We can now plot the expression of some marker genes to see their
distribution in all the clusters

``` r
VlnPlot(pbmc10_05, features = c("Fabp4", "Gsn"))
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-112-1.png)<!-- -->

From the plot we can see that Fabp4 is expressed at high level only in
cells of the first cluster, while Gsn is expressed at highest level only
in the cells of the cluster 1 and with a high value in some of the cells
of the cluster 2.

``` r
DotPlot(pbmc10_05, features = c("Ifit1", "Lyz2", "Acta1", "Smoc2", "Cd79a", "Fmod", "Rgs5", "Sdc4", "Cd3g"))
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-113-1.png)<!-- -->

From the dotplot is visible that each of the plotted genes is a unique
marker of a cluster: in fact even if Sdc4 is expressed also in cells of
other clusters, each marker gene is expressed with the highest average
value and in the highest percentage of cells only of one cluster.

We can plot the UMAP for each marker gene where the cells expressing
that gene are colored

``` r
FeaturePlot(pbmc10_05, features = c("Ifit1", "Lyz2", "Acta1", "Smoc2", "Cd79a", "Fmod", "Rgs5", "Sdc4", "Cd3g"))
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-114-1.png)<!-- -->

As before, we can see that quite all the markers are expressed only in
one cluster of cells with the exception of Ifit1 and Sdc4 that - as
already seen in the dotplot - even if expressed with highest value by
all the cells of one cluster are also expressed by a smaller percentage
of cells of other clusters.

Eventually we can **plot the cells with the corresponding cell type**

``` r
new.cluster.ids <- c("endothelial cell", "mesenchymal stem cell", "fibroblast", "skeletal muscle satellite cell", "B", "macrophage", "T", "smooth muscle cell", "myocytes")
names(new.cluster.ids) <- levels(pbmc10_05)
pbmc10_05 <- RenameIdents(pbmc10_05, new.cluster.ids)
DimPlot(pbmc10_05, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-115-1.png)<!-- -->

To conclude we can see that the result of this analysis - performed by
filtering more cells, finding in a different way the most variable genes
and considering the first 10 principal components and a resolution
parameter of 0.5 - is very different from the previous ones. First of
all this analysis returns only 9 cell types - Schwann cells are not
identified - and what is very strange is the fact that the cluster of
endothelial cells is made of a big concentrated cluster plus a small
cluster far from it and very near to mesenchymal stem cells and smooth
muscle cells.

We also plot this result using tsne visualization

``` r
new.cluster.ids <- c("endothelial cell", "mesenchymal stem cell", "fibroblast", "skeletal muscle satellite cell", "B", "macrophage", "T", "smooth muscle cell", "myocytes")
names(new.cluster.ids) <- levels(pbmc10_05)
pbmc10_05 <- RenameIdents(pbmc10_05, new.cluster.ids)
DimPlot(pbmc10_05, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
```

![](final-scRNA-seq_files/figure-gfm/unnamed-chunk-116-1.png)<!-- -->

Our result obtained using the article parameters is quite different -
both because the position and size of the different clusters is
different from the article’s one and also because we obtained
fibroblast, myocytes and smooth muscle cell types that are not present
in the article’s result. This can be explained by the fact that **the
data that the authors of the article used is different from the data on
which we performed the analysis: in particular their data contains less
cells and a lower overall number of features.**

Anyway we decided to perform and show the analysis using the parameters
and modifications suggested by the article, even if our dataset is
different, because we wanted to see whether our data was not too
different and so if we could achieve the same result. Unfortunately this
is not the case!

Eventually we can also see that this result is very different from the
one obtained looking at the first 16 PCs and using a resolution of 0.2
when performing the clustering. In particular this analysis is less
precise since it was not able to detect Schwann cell type and
furthermore we can see from the UMAP that the clusters are less well
separated than the result obtained using 16 PCs and a resolution of 0.2.

In conclusion **the best parameters - that allow to detect more cell
types from the muscle sample associated with more defined clusters - are
PC=16 and a resolution of 0.2.**
