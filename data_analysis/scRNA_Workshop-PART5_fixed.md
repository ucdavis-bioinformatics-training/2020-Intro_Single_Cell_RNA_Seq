---
title: "Single Cell RNAseq Part 5"
author: "Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

## Load libraries

```r
library(Seurat)
library(ggplot2)
```

## Load the Seurat object

```r
load(file="pca_sample_corrected.RData")
experiment.aggregate
```

<div class='r_output'> An object of class Seurat 
 12811 features across 2681 samples within 1 assay 
 Active assay: RNA (12811 features, 2000 variable features)
  1 dimensional reduction calculated: pca
</div>
## Identifying clusters

Seurat implements an graph-based clustering approach. Distances between the cells are calculated based on previously identified PCs. Seurat approach was heavily inspired by recent manuscripts which applied graph-based clustering approaches to scRNAseq data. Briefly, Seurat identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. First calculate k-nearest neighbors (KNN) and construct the SNN graph. Then optimize the modularity function to determine clusters. For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B.

The FindClusters function implements the procedure, and contains a resolution parameter that sets the granularity of the downstream clustering, with increased values leading to a greater number of clusters. I tend to like to perform a series of resolutions, investigate and choose.


```r
use.pcs = 1:29 

?FindNeighbors
experiment.aggregate <- FindNeighbors(experiment.aggregate, reduction="pca", dims = use.pcs)
```

<div class='r_output'> Computing nearest neighbor graph
</div>
<div class='r_output'> Computing SNN
</div>
```r
?FindCluster
```

<div class='r_output'> No documentation for 'FindCluster' in specified packages and libraries:
 you could try '??FindCluster'
</div>
```r
experiment.aggregate <- FindClusters(
    object = experiment.aggregate, 
    resolution = seq(0.25,4,0.25), 
    verbose = FALSE
)
```

Lets first investigate how many clusters each resolution produces and set it to the smallest resolutions of 0.5 (fewest clusters). 


```r
sapply(grep("res",colnames(experiment.aggregate@meta.data),value = TRUE),
       function(x) length(unique(experiment.aggregate@meta.data[,x])))
```

<div class='r_output'> RNA_snn_res.0.25  RNA_snn_res.0.5 RNA_snn_res.0.75    RNA_snn_res.1 
               10               13               15               17 
 RNA_snn_res.1.25  RNA_snn_res.1.5 RNA_snn_res.1.75    RNA_snn_res.2 
               18               20               24               24 
 RNA_snn_res.2.25  RNA_snn_res.2.5 RNA_snn_res.2.75    RNA_snn_res.3 
               25               26               27               28 
 RNA_snn_res.3.25  RNA_snn_res.3.5 RNA_snn_res.3.75    RNA_snn_res.4 
               28               28               28               30
</div>
```r
Idents(experiment.aggregate) <- "RNA_snn_res.0.5"
```

Finally,  lets produce a table of cluster to sample assignments.

```r
table(Idents(experiment.aggregate),experiment.aggregate$orig.ident)
```

<div class='r_output'>     
      UCD_Adj_VitE UCD_Supp_VitE UCD_VitE_Def
   0           156           171          154
   1           119           140          107
   2            78           118          160
   3            93            81           93
   4            71            82           66
   5            55            77           86
   6            60            67           65
   7            50            80           62
   8            23            45           39
   9            34            33           38
   10           31            19           18
   11           18            20           19
   12           20            14           19
</div>
tSNE dimensionality reduction plots are then used to visualise clustering results. As input to the tSNE, you should use the same PCs as input to the clustering analysis.


```r
experiment.aggregate <- RunTSNE(
  object = experiment.aggregate,
  reduction.use = "pca",
  dims.use = use.pcs,
  do.fast = TRUE)
```

Plot TSNE coloring by the slot 'ident' (default).

```r
DimPlot(object = experiment.aggregate, pt.size=0.5, reduction = "tsne", label = T)
```

<div class='r_output'> Warning: Using `as.character()` on a quosure is deprecated as of rlang 0.3.0.
 Please use `as_label()` or `as_name()` instead.
 This warning is displayed once per session.
</div>
![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-7-1.png)<!-- -->


Plot TSNE coloring by the clustering resolution 4

```r
DimPlot(object = experiment.aggregate, group.by="RNA_snn_res.4", pt.size=0.5, reduction = "tsne", label = T)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

uMAP dimensionality reduction plot.


```r
experiment.aggregate <- RunUMAP(
  object = experiment.aggregate,
  dims = use.pcs)
```

<div class='r_output'> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
 To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
 This message will be shown once per session
</div>
<div class='r_output'> 06:42:45 UMAP embedding parameters a = 0.9922 b = 1.112
</div>
<div class='r_output'> 06:42:45 Read 2681 rows and found 29 numeric columns
</div>
<div class='r_output'> 06:42:45 Using Annoy for neighbor search, n_neighbors = 30
</div>
<div class='r_output'> 06:42:45 Building Annoy index with metric = cosine, n_trees = 50
</div>
<div class='r_output'> 0%   10   20   30   40   50   60   70   80   90   100%
</div>
<div class='r_output'> [----|----|----|----|----|----|----|----|----|----|
</div>
<div class='r_output'> **************************************************|
 06:42:46 Writing NN index file to temp file /var/folders/74/h45z17f14l9g34tmffgq9nkw0000gn/T//Rtmp6y8x3z/file586337bc35b
 06:42:46 Searching Annoy index using 1 thread, search_k = 3000
 06:42:46 Annoy recall = 100%
 06:42:47 Commencing smooth kNN distance calibration using 1 thread
 06:42:47 Initializing from normalized Laplacian + noise
 06:42:47 Commencing optimization for 500 epochs, with 107472 positive edges
 06:42:50 Optimization finished
</div>
Plot uMap coloring by the slot 'ident' (default).

```r
DimPlot(object = experiment.aggregate, pt.size=0.5, reduction = "umap", label = T)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


FeaturePlot can be used to color cells with a 'feature', non categorical data, like number of UMIs

```r
FeaturePlot(experiment.aggregate, features = c('nCount_RNA'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-11-1.png)<!-- -->
and number of genes present

```r
FeaturePlot(experiment.aggregate, features = c('nFeature_RNA'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

percent mitochondrial 

```r
FeaturePlot(experiment.aggregate, features = c('percent.mito'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

TSNE plot by cell cycle

```r
DimPlot(object = experiment.aggregate, pt.size=0.5, group.by = "cell.cycle", reduction = "tsne" )
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-14-1.png)<!-- -->


## Building  a  tree relating the 'average' cell from each cluster. Tree is estimated based on a distance matrix constructed in either gene expression space or PCA space.


```r
experiment.aggregate <- BuildClusterTree(
  experiment.aggregate, dims = use.pcs)

PlotClusterTree(experiment.aggregate)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-15-1.png)<!-- -->


```r
DimPlot(object = experiment.aggregate, pt.size=0.5, label = TRUE, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

Merge Clustering results

```r
experiment.merged = experiment.aggregate
# originally set clusters to resolutionm 0.5
Idents(experiment.merged) <- "RNA_snn_res.0.5"

table(Idents(experiment.merged))
```

<div class='r_output'> 
   0   1   2   3   4   5   6   7   8   9  10  11  12 
 481 366 356 267 219 218 192 192 107 105  68  57  53
</div>
```r
# based on TSNE and Heirarchical tree
# merge clusters 6 and 7 into 0 and cluster 9 into 13
experiment.merged <- RenameIdents(
  object = experiment.merged,
  '6' = '0', '7' = '0', '9' = '13'
)

table(Idents(experiment.merged))
```

<div class='r_output'> 
   0  13   1   2   3   4   5   8  10  11  12 
 865 105 366 356 267 219 218 107  68  57  53
</div>
```r
DimPlot(object = experiment.merged, pt.size=0.5, label = T, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

```r
experiment.examples <- experiment.merged
# in order to reporder the clusters for plotting purposes
# take a look at the levels, which indicates the ordering
levels(experiment.examples@active.ident)
```

<div class='r_output'>  [1] "0"  "13" "1"  "2"  "3"  "4"  "5"  "8"  "10" "11" "12"
</div>
```r
# relevel setting 5 to the first factor
experiment.examples@active.ident <- relevel(experiment.examples@active.ident, "5")
levels(experiment.examples@active.ident)
```

<div class='r_output'>  [1] "5"  "0"  "13" "1"  "2"  "3"  "4"  "8"  "10" "11" "12"
</div>
```r
# now cluster 5 is the "first" factor

# relevel all the factors to the order I want
Idents(experiment.examples) <- factor(experiment.examples@active.ident, levels=c("5","13","1","2","3","0","4","8","11","12","10","14"))
levels(experiment.examples@active.ident)
```

<div class='r_output'>  [1] "5"  "13" "1"  "2"  "3"  "0"  "4"  "8"  "11" "12" "10"
</div>
```r
DimPlot(object = experiment.examples, pt.size=0.5, label = T, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-17-2.png)<!-- -->

```r
### Re-assign clustering result to resolution 4 for cells in cluster 0 (@ reslution 0.5) [adding a R prefix]
newIdent = as.character(Idents(experiment.examples))
newIdent[newIdent == '0'] = paste0("R",as.character(experiment.examples$RNA_snn_res.4[newIdent == '0']))

Idents(experiment.examples) <- as.factor(newIdent)

table(Idents(experiment.examples))
```

<div class='r_output'> 
   1  10  11  12  13   2   3   4   5   8  R0 R10 R11 R12 R14 R15 R17 R20 R22 R26 
 366  68  57  53 105 356 267 219 218 107 176  20   2  89   1  74   3  61   1   1 
 R27 R28 R29  R5  R7  R8  R9 
  41  33   1 143 118   1 100
</div>
```r
DimPlot(object = experiment.examples, pt.size=0.5, label = T, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-17-3.png)<!-- -->

Plot TSNE coloring by the slot 'orig.ident' (sample names) with alpha colors turned on.

```r
DimPlot(object = experiment.aggregate, group.by="orig.ident", pt.size=0.5, reduction = "tsne" )
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

```r
## Pretty tsne using alpha
p <- DimPlot(object = experiment.aggregate, group.by="orig.ident", pt.size=0.5, reduction = "tsne")
alpha.use <- 2/5
p$layers[[1]]$mapping$alpha <- alpha.use
p + scale_alpha_continuous(range = alpha.use, guide = F)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-18-2.png)<!-- -->

Removing cells assigned to clusters from a plot, So here plot all clusters but clusters "3" and "5"

```r
# create a new tmp object with those removed 
experiment.aggregate.tmp <- experiment.aggregate[,-which(Idents(experiment.aggregate) %in% c("3","5"))]

dim(experiment.aggregate)
```

<div class='r_output'> [1] 12811  2681
</div>
```r
dim(experiment.aggregate.tmp)
```

<div class='r_output'> [1] 12811  2196
</div>
```r
DimPlot(object = experiment.aggregate.tmp, group.by="orig.ident", pt.size=0.5, reduction = "tsne", label = T)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

## Identifying Marker Genes

Seurat can help you find markers that define clusters via differential expression.

`FindMarkers` identifies markers for a cluster relative to all other clusters.

`FindAllMarkers` does so for all clusters

`FindAllMarkersNode` defines all markers that split a Node __(Warning: need to validate)__


```r
?FindMarkers

markers = FindMarkers(experiment.merged, ident.1=c(10), genes.use = VariableFeatures(experiment.merged))

head(markers)
```

<div class='r_output'>                 p_val avg_logFC pct.1 pct.2     p_val_adj
 Foxp2   3.059052e-146 1.1096456 0.324 0.002 3.918952e-142
 Abi3bp   9.558391e-66 0.9537030 0.206 0.004  1.224525e-61
 Calb1    2.923569e-54 1.6284265 0.294 0.015  3.745385e-50
 Slc10a4  6.562276e-47 0.3783238 0.103 0.001  8.406931e-43
 Kcnk9    1.347345e-45 0.5768638 0.118 0.002  1.726084e-41
 Col4a4   5.955819e-35 0.5071647 0.103 0.002  7.630000e-31
</div>
```r
dim(markers)
```

<div class='r_output'> [1] 1299    5
</div>
```r
table(markers$avg_logFC > 0)
```

<div class='r_output'> 
 FALSE  TRUE 
   735   564
</div>
 
pct.1 and pct.2 are the proportion of cells with expression above 0 in ident.1 and ident.2 respectively. p_val is the raw p_value associated with the differntial expression test with adjusted value in p_val_adj. avg_logFC is the average log fold change difference between the two groups. 
 
avg_diff (lines 130, 193 and) appears to be the difference in log(x = mean(x = exp(x = x) - 1) + 1) between groups.  It doesn’t seem like this should work out to be the signed ratio of pct.1 to pct.2 so I must be missing something.  It doesn’t seem to be related at all to how the p-values are calculated so maybe it doesn’t matter so much, and the sign is probably going to be pretty robust to how expression is measured.

Can use a violin plot to visualize the expression pattern of some markers

```r
VlnPlot(object = experiment.merged, features = rownames(markers)[1:2], pt.size = 0.05)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

Or a feature plot

```r
FeaturePlot(
    experiment.merged, 
    head(rownames(markers), n=6), 
    cols = c("lightgrey", "blue"), 
    ncol = 2
)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

```r
FeaturePlot(    
    experiment.merged, 
    "Fxyd1", 
    cols = c("lightgrey", "blue") 
)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-22-2.png)<!-- -->

FindAllMarkers can be used to automate the process across all genes.
__WARNING: TAKES A LONG TIME TO RUN__


```r
markers_all <- FindAllMarkers(
    object = experiment.merged, 
    only.pos = TRUE, 
    min.pct = 0.25, 
    thresh.use = 0.25
)
```

<div class='r_output'> Calculating cluster 0
</div>
<div class='r_output'> Calculating cluster 13
</div>
<div class='r_output'> Calculating cluster 1
</div>
<div class='r_output'> Calculating cluster 2
</div>
<div class='r_output'> Calculating cluster 3
</div>
<div class='r_output'> Calculating cluster 4
</div>
<div class='r_output'> Calculating cluster 5
</div>
<div class='r_output'> Calculating cluster 8
</div>
<div class='r_output'> Calculating cluster 10
</div>
<div class='r_output'> Calculating cluster 11
</div>
<div class='r_output'> Calculating cluster 12
</div>
```r
dim(markers_all)
```

<div class='r_output'> [1] 3896    7
</div>
```r
head(markers_all)
```

<div class='r_output'>                       p_val avg_logFC pct.1 pct.2     p_val_adj cluster
 Fxyd7         1.506978e-112 1.3535140 0.521 0.131 1.930589e-108       0
 Marcks         1.381731e-76 0.8007480 0.721 0.372  1.770135e-72       0
 Ppp3ca         3.782331e-69 0.6427877 0.928 0.826  4.845544e-65       0
 Pcp4           1.610965e-67 0.8076783 0.554 0.225  2.063807e-63       0
 6330403K07Rik  3.056161e-57 0.7862346 0.592 0.305  3.915248e-53       0
 Tagln3         8.962286e-54 0.7980438 0.551 0.286  1.148158e-49       0
                        gene
 Fxyd7                 Fxyd7
 Marcks               Marcks
 Ppp3ca               Ppp3ca
 Pcp4                   Pcp4
 6330403K07Rik 6330403K07Rik
 Tagln3               Tagln3
</div>
```r
table(table(markers_all$gene))
```

<div class='r_output'> 
    1    2    3    4    5 
 1723  594  232   66    5
</div>
```r
markers_all_single <- markers_all[markers_all$gene %in% names(table(markers_all$gene))[table(markers_all$gene) == 1],]

dim(markers_all_single)
```

<div class='r_output'> [1] 1723    7
</div>
```r
table(table(markers_all_single$gene))
```

<div class='r_output'> 
    1 
 1723
</div>
```r
table(markers_all_single$cluster)
```

<div class='r_output'> 
   0  13   1   2   3   4   5   8  10  11  12 
  24  70 116 232 162 690 141  57  15 123  93
</div>
```r
head(markers_all_single)
```

<div class='r_output'>                 p_val avg_logFC pct.1 pct.2    p_val_adj cluster     gene
 S100a13  2.655811e-50 0.6909175 0.655 0.436 3.402360e-46       0  S100a13
 Gm13889  1.233554e-44 0.6728831 0.683 0.475 1.580305e-40       0  Gm13889
 Nrsn1    5.173753e-44 0.7934783 0.703 0.487 6.628095e-40       0    Nrsn1
 Map1lc3a 4.422443e-27 0.2803135 0.924 0.858 5.665592e-23       0 Map1lc3a
 Ywhaq    1.120435e-22 0.3145009 0.827 0.759 1.435390e-18       0    Ywhaq
 Nrip1    3.862595e-20 0.6182926 0.332 0.195 4.948370e-16       0    Nrip1
</div>
Plot a heatmap of genes by cluster for the top 5 marker genes per cluster

```r
library(dplyr)
```

<div class='r_output'> 
 Attaching package: 'dplyr'
</div>
<div class='r_output'> The following objects are masked from 'package:stats':
 
     filter, lag
</div>
<div class='r_output'> The following objects are masked from 'package:base':
 
     intersect, setdiff, setequal, union
</div>
```r
top5 <- markers_all_single %>% group_by(cluster) %>% top_n(5, avg_logFC)
dim(top5)
```

<div class='r_output'> [1] 55  7
</div>
```r
DoHeatmap(
    object = experiment.merged, 
    features = top5$gene
) 
```

<div class='r_output'> Warning in DoHeatmap(object = experiment.merged, features = top5$gene): The
 following features were omitted as they were not found in the scale.data slot
 for the RNA assay: Mest, Nwd2, Ctnnd2, Nrip1, Nrsn1, Gm13889
</div>
![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-24-1.png)<!-- -->


```r
# Get expression of genes for cells in and out of each cluster
getGeneClusterMeans <- function(gene, cluster){
  x <- GetAssayData(experiment.merged)[gene,]
  m <- tapply(x, ifelse(Idents(experiment.merged) == cluster, 1, 0), mean)
  mean.in.cluster <- m[2]
  mean.out.of.cluster <- m[1]
  return(list(mean.in.cluster = mean.in.cluster, mean.out.of.cluster = mean.out.of.cluster))
}

## for sake of time only using first six (head)
means <- mapply(getGeneClusterMeans, head(markers_all[,"gene"]), head(markers_all[,"cluster"]))
means <- matrix(unlist(means), ncol = 2, byrow = T)

colnames(means) <- c("mean.in.cluster", "mean.out.of.cluster")
rownames(means) <- head(markers_all[,"gene"])
markers_all2 <- cbind(head(markers_all), means)
head(markers_all2)
```

<div class='r_output'>                       p_val avg_logFC pct.1 pct.2     p_val_adj cluster
 Fxyd7         1.506978e-112 1.3535140 0.521 0.131 1.930589e-108       0
 Marcks         1.381731e-76 0.8007480 0.721 0.372  1.770135e-72       0
 Ppp3ca         3.782331e-69 0.6427877 0.928 0.826  4.845544e-65       0
 Pcp4           1.610965e-67 0.8076783 0.554 0.225  2.063807e-63       0
 6330403K07Rik  3.056161e-57 0.7862346 0.592 0.305  3.915248e-53       0
 Tagln3         8.962286e-54 0.7980438 0.551 0.286  1.148158e-49       0
                        gene mean.in.cluster mean.out.of.cluster
 Fxyd7                 Fxyd7        1.402763           0.2900659
 Marcks               Marcks        1.604625           0.7205869
 Ppp3ca               Ppp3ca        2.748383           2.0118646
 Pcp4                   Pcp4        1.470163           0.5574350
 6330403K07Rik 6330403K07Rik        1.354101           0.5966355
 Tagln3               Tagln3        1.210059           0.5233606
</div>
## Finishing up clusters.

At this point in time you should use the tree, markers, domain knowledge, and goals to finalize your clusters. This may mean adjusting PCA to use, mergers clusters together, choosing a new resolutions, etc. When finished you can further name it cluster by something more informative. Ex.

```r
experiment.clusters <- experiment.aggregate
experiment.clusters <- RenameIdents(
  object = experiment.clusters,
  '0' = 'cell_type_A',
  '1' = 'cell_type_B',
  '2' = 'cell_type_C'
)
# and so on

DimPlot(object = experiment.clusters, pt.size=0.5, label = T, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-26-1.png)<!-- -->

```r
experiment.merged$finalcluster <- Idents(experiment.merged)
```

## Subsetting samples
If you want to look at the representation of just one sample, or sets of samples

```r
experiment.sample2 <- subset(experiment.merged, orig.ident == "UCD_Supp_VitE")

DimPlot(object = experiment.sample2, group.by = "RNA_snn_res.0.5", pt.size=0.5, label = TRUE, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

```r
FeaturePlot(experiment.sample2, features =c('Calca'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-27-2.png)<!-- -->

```r
FeaturePlot(experiment.sample2, features =c('Adcyap1'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-27-3.png)<!-- -->

```r
experiment.batch1 <- subset(experiment.merged, batchid == "Batch1")

DimPlot(object = experiment.batch1, group.by = "RNA_snn_res.0.5", pt.size=0.5, label = TRUE, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-27-4.png)<!-- -->

### Adding in a new metadata column representing samples within clusters


```r
experiment.merged$samplecluster = paste(experiment.merged$orig.ident,experiment.merged$finalcluster,sep = '-')

# set the identity to the new variable 
Idents(experiment.merged) <- "samplecluster"

markers.comp <- FindMarkers(experiment.merged, ident.1 = "UCD_Adj_VitE-0", ident.2= c("UCD_Supp_VitE-0","UCD_VitE_Def-0"))

markers.comp
```

<div class='r_output'>                 p_val  avg_logFC pct.1 pct.2    p_val_adj
 Actb     8.844395e-11 -0.5646784 0.996 0.982 1.133055e-06
 Tmsb10   3.412356e-10  0.2598487 0.989 0.950 4.371569e-06
 Rpl21    3.055635e-09  0.3875170 0.902 0.735 3.914574e-05
 Rpl23a   1.736979e-07  0.2774184 0.936 0.813 2.225243e-03
 Pcp4     2.426460e-07  0.4390500 0.684 0.496 3.108538e-03
 Rpl17    1.108312e-06  0.2637329 0.883 0.751 1.419859e-02
 Snrpg    7.328599e-06  0.2521357 0.259 0.132 9.388668e-02
 S100a11  1.767595e-05  0.2932752 0.286 0.157 2.264466e-01
 BC005624 3.637158e-05  0.3149532 0.282 0.165 4.659563e-01
 Atp5e    4.679088e-05  0.2663626 0.718 0.618 5.994380e-01
 Wdfy1    4.759243e-05 -0.2946807 0.060 0.160 6.097066e-01
 Gpx3     8.375777e-05  0.3202155 0.357 0.222 1.000000e+00
 Arhgap15 9.299868e-05  0.2873916 0.233 0.127 1.000000e+00
 Pdap1    9.567528e-05  0.2508146 0.447 0.299 1.000000e+00
 Dbpht2   5.669139e-04  0.2673831 0.316 0.205 1.000000e+00
 Mt2      7.474565e-04  0.3494247 0.188 0.104 1.000000e+00
 Cbx3     1.442723e-03 -0.3846061 0.342 0.431 1.000000e+00
 Cartpt   2.145946e-03  0.3111675 0.184 0.105 1.000000e+00
 Vdac1    3.575640e-03 -0.3064423 0.526 0.568 1.000000e+00
 Fam168b  4.024716e-03 -0.2918939 0.165 0.240 1.000000e+00
 Abhd2    1.112627e-02 -0.5929060 0.402 0.454 1.000000e+00
 Eno2     1.556555e-02 -0.3043541 0.421 0.462 1.000000e+00
 Tagln2   1.986845e-02 -0.2745054 0.098 0.154 1.000000e+00
 Dcaf12   2.285967e-02 -0.2503752 0.150 0.210 1.000000e+00
 Zbtb4    2.408377e-02 -0.2736469 0.248 0.302 1.000000e+00
 Arf3     2.703628e-02 -0.2759801 0.571 0.574 1.000000e+00
 Usp22    2.758262e-02 -0.2705576 0.305 0.362 1.000000e+00
 Itm2c    3.007439e-02 -0.2978899 0.308 0.354 1.000000e+00
 S100a16  3.224776e-02 -0.2978945 0.098 0.149 1.000000e+00
 Nacc2    3.357621e-02 -0.2599636 0.477 0.504 1.000000e+00
 Necab1   5.669348e-02 -0.2986063 0.158 0.207 1.000000e+00
 Birc6    6.090972e-02 -0.2689956 0.158 0.204 1.000000e+00
 Nefh     6.359203e-02 -0.2753826 0.256 0.306 1.000000e+00
 S100b    7.647967e-02 -0.2557280 0.214 0.265 1.000000e+00
 Etv5     7.792293e-02 -0.2811197 0.184 0.225 1.000000e+00
 Clasp2   9.020684e-02 -0.2664070 0.207 0.252 1.000000e+00
 Hdlbp    9.154369e-02 -0.2878492 0.177 0.215 1.000000e+00
 Nfia     1.008055e-01 -0.3377224 0.203 0.242 1.000000e+00
 AI593442 1.025835e-01 -0.2862616 0.252 0.289 1.000000e+00
 Mt1      1.745423e-01  0.3141396 0.485 0.441 1.000000e+00
 Pam      1.851790e-01 -0.3745617 0.511 0.526 1.000000e+00
 Ythdf2   1.890034e-01 -0.2655538 0.342 0.354 1.000000e+00
 Nat8l    1.969270e-01 -0.2708323 0.259 0.284 1.000000e+00
 Cntn1    1.970727e-01 -0.3056274 0.308 0.334 1.000000e+00
 Alkbh5   2.100413e-01 -0.2519123 0.211 0.242 1.000000e+00
 Ehd3     2.123385e-01 -0.2710123 0.192 0.215 1.000000e+00
 Ddhd1    2.200160e-01 -0.2688655 0.195 0.220 1.000000e+00
 Spry2    2.782289e-01 -0.3101408 0.211 0.225 1.000000e+00
 Bhlhe41  3.928239e-01 -0.2563678 0.365 0.374 1.000000e+00
 Wtap     4.174552e-01 -0.2745871 0.372 0.374 1.000000e+00
</div>
```r
experiment.subset <- subset(experiment.merged, samplecluster %in%  c( "UCD_Adj_VitE-0", "UCD_Supp_VitE-0" ))
DoHeatmap(experiment.subset, features = rownames(markers.comp))
```

<div class='r_output'> Warning in DoHeatmap(experiment.subset, features = rownames(markers.comp)): The
 following features were omitted as they were not found in the scale.data slot
 for the RNA assay: Wtap, Ddhd1, Ythdf2, Hdlbp, Clasp2, Birc6, Necab1, Itm2c,
 Usp22, Arf3, Zbtb4, Dcaf12, Fam168b, Vdac1, Pdap1, Arhgap15, Wdfy1, Atp5e,
 BC005624, Snrpg, Rpl17, Rpl23a, Rpl21, Tmsb10
</div>
![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

```r
Idents(experiment.merged) <- "finalcluster"
```

And last lets save all the objects in our session.

```r
save(list=ls(), file="clusters_seurat_object.RData")
```

## Get the next Rmd file

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/master/scrnaseq_analysis/scRNA_Workshop-PART6.Rmd", "scRNA_Workshop-PART6.Rmd")
```

## Session Information

```r
sessionInfo()
```

<div class='r_output'> R version 4.0.0 (2020-04-24)
 Platform: x86_64-apple-darwin17.0 (64-bit)
 Running under: macOS Catalina 10.15.4
 
 Matrix products: default
 BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
 LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
 
 locale:
 [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
 
 attached base packages:
 [1] stats     graphics  grDevices datasets  utils     methods   base     
 
 other attached packages:
 [1] dplyr_0.8.5   ggplot2_3.3.0 Seurat_3.1.5 
 
 loaded via a namespace (and not attached):
  [1] httr_1.4.1          tidyr_1.0.3         jsonlite_1.6.1     
  [4] viridisLite_0.3.0   splines_4.0.0       leiden_0.3.3       
  [7] assertthat_0.2.1    BiocManager_1.30.10 renv_0.10.0        
 [10] yaml_2.2.1          ggrepel_0.8.2       globals_0.12.5     
 [13] pillar_1.4.4        lattice_0.20-41     limma_3.44.1       
 [16] glue_1.4.1          reticulate_1.15     digest_0.6.25      
 [19] RColorBrewer_1.1-2  colorspace_1.4-1    cowplot_1.0.0      
 [22] htmltools_0.4.0     Matrix_1.2-18       plyr_1.8.6         
 [25] pkgconfig_2.0.3     tsne_0.1-3          listenv_0.8.0      
 [28] purrr_0.3.4         patchwork_1.0.0     scales_1.1.1       
 [31] RANN_2.6.1          RSpectra_0.16-0     Rtsne_0.15         
 [34] tibble_3.0.1        farver_2.0.3        ellipsis_0.3.1     
 [37] withr_2.2.0         ROCR_1.0-11         pbapply_1.4-2      
 [40] lazyeval_0.2.2      survival_3.1-12     magrittr_1.5       
 [43] crayon_1.3.4        evaluate_0.14       future_1.17.0      
 [46] nlme_3.1-147        MASS_7.3-51.5       ica_1.0-2          
 [49] tools_4.0.0         fitdistrplus_1.1-1  data.table_1.12.8  
 [52] lifecycle_0.2.0     stringr_1.4.0       plotly_4.9.2.1     
 [55] munsell_0.5.0       cluster_2.1.0       irlba_2.3.3        
 [58] compiler_4.0.0      rsvd_1.0.3          rlang_0.4.6        
 [61] grid_4.0.0          ggridges_0.5.2      RcppAnnoy_0.0.16   
 [64] htmlwidgets_1.5.1   igraph_1.2.5        labeling_0.3       
 [67] rmarkdown_2.1       gtable_0.3.0        codetools_0.2-16   
 [70] reshape2_1.4.4      R6_2.4.1            gridExtra_2.3      
 [73] zoo_1.8-8           knitr_1.28          uwot_0.1.8         
 [76] future.apply_1.5.0  KernSmooth_2.23-16  ape_5.3            
 [79] stringi_1.4.6       parallel_4.0.0      Rcpp_1.0.4.6       
 [82] vctrs_0.3.0         sctransform_0.2.1   png_0.1-7          
 [85] tidyselect_1.1.0    xfun_0.13           lmtest_0.9-37
</div>