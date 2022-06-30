# scRNA-seq-analysis

The analysis is performed using a set of scRNA-seq data - produced by 10x Genomics Chromium - dowloaded from PanglaoDB (https://panglaodb.se/view_data.php?sra=SRA653146&srs=SRS3044259).
These data derive from 2737 cells of a mouse muscle sample.

The analysis aims at identifying which cell types are present in the data using different combinations of parameters to perform the clustering - 
the number of Principal Components used to extract only the the most variable genes expressed by the cells and resolution parameter that controls the granularity
of the clusters.
Furthermore the analysis is repeated using two different ways of selecting the most variable genes.
