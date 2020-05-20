---
title: "Single Cell RNAseq Part 1"
author: "Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

# Single Cell Analysis with Seurat and some custom code!

[Seurat](http://satijalab.org/seurat/) is a popular R package that is designed for QC, analysis, and exploration of single cell data. Seurat aims to enable users to identify and interpret sources of heterogeneity from single cell transcriptomic measurements, and to integrate diverse types of single cell data. Further, the authors provide several [tutorials](https://satijalab.org/seurat/vignettes.html), on their website.

The expression_tables_cellrangerV3.zip file contains the single cell matrix files for the three samples. These are isolated mouse cells ran on the 10X genomics platform (3' expression V2) for single cell RNA sequencing, sequenced with UC Davis on 1 HiSeq 4000 lane. The Experiment contains 3 samples, each merged from 2 original samples and then randomly subsamples to 1000 cells each.

The three samples are, Dorsal root ganglion neurons :
At weaning, Ttpa+/+ mice were fed a normal diet (35 mg of dl-α-tocopheryl acetate/kg, vitE+) while and Ttpa-/- mice were fed either an α-TOH-deficient diet (<10 mg of dl-α-tocopheryl acetate/kg, vitE-), or α-TOH-supplemented diet (600 mg of dl-α-tocopheryl acetate/kg, vitE+++) diet.

* UCD_Adj_VitE - normal + Vitamin E
* UCD_Supp_VitE - Vitamin E supplimented by diet.
* UCD_VitE_Def - Vitamin E deficient animals

We start with loading needed libraries for R, at this time all we need is the package [Seurat](http://satijalab.org/seurat/).

```r
library(Seurat)
library(biomaRt)
```

## Load the Cell Ranger Matrix Data and create the base Seurat object.
Cell Ranger provides a function `cellranger aggr` that will combine multiple samples into a single matrix file. However, when processing data in R and Seurat this is unnecessary and we can aggregate them in R.

Seurat provides a function `Read10X` to read in 10X data folder. First we read in data from each individual sample folder. First, we initialize the Seurat object (`CreateSeuratObject`) with the raw (non-normalized data). Keep all genes expressed in >= 10 cells. Keep all cells with at least 200 detected genes. Also extracting sample names, calculating and adding in the metadata mitochondrial percentage of each cell. Adding in the metadata batchid and cell cycle. Finally, saving the raw Seurat object.


```r
dataset_loc <- "./expression_tables_cellrangerV3"
ids <- c("UCD_Adj_VitE", "UCD_Supp_VitE", "UCD_VitE_Def")

d10x.data <- sapply(ids, function(i){
  d10x <- Read10X(file.path(dataset_loc,i,"outs/filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})

experiment.data <- do.call("cbind", d10x.data)

experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = "scRNA workshop",
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "\\-")
```

### The percentage of reads that map to the mitochondrial genome

* Low-quality / dying cells often exhibit extensive mitochondrial contamination.
* We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of features.
* We use the set of all genes, in mouse these genes can be identified as those that begin with 'mt', in human data they begin with MT.


```r
experiment.aggregate$percent.mito <- PercentageFeatureSet(experiment.aggregate, pattern = "^mt-")
```


## Calculate cell cycle, add to meta data
Using the package scran, get the mouse cell cycle markers and a mapping of m
**There may be issues with timeine for the Biomart server, just keep trying**

```r
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
# Convert to matrix for use in cycle
mat <- as.matrix(GetAssayData(experiment.aggregate))

# Convert rownames to ENSEMBL IDs, Using biomaRt
ensembl<- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
```

<div class='r_output'> Ensembl site unresponsive, trying uswest mirror
</div>
```r
anno <- getBM( values=rownames(mat), attributes=c("mgi_symbol","ensembl_gene_id") , filters= "mgi_symbol"  ,mart=ensembl)

ord <- match(rownames(mat), anno$mgi_symbol) # use anno$mgi_symbol if via biomaRt
rownames(mat) <- anno$ensembl_gene_id[ord] # use anno$ensembl_gene_id if via biomaRt
drop <- which(is.na(rownames(mat)))
mat <- mat[-drop,]
cycles <- scran::cyclone(mat, pairs=mm.pairs)
tmp <- data.frame(cell.cycle = cycles$phases)
rownames(tmp) <- colnames(mat)
experiment.aggregate <- AddMetaData(experiment.aggregate, tmp)
```

### Cell-Cycle with Seurat, the list of genes comes with Seurat (only for human)

First need to convert to mouse symbols, we'll use Biomart for that too.

```r
table(Idents(experiment.aggregate))
```

<div class='r_output'> 
  UCD_Adj_VitE UCD_Supp_VitE  UCD_VitE_Def 
           904          1000           992
</div>
```r
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
   
  genes = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
   
  humanx <- unique(genes[, 2])
   
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
 
m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
```

<div class='r_output'> [1] "Gmnn"  "Rad51" "Cdca7" "Prim1" "Slbp"  "Mcm7"
</div>
```r
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)
```

<div class='r_output'> Ensembl site unresponsive, trying uswest mirror
</div>
<div class='r_output'> [1] "Tacc3"  "Tmpo"   "Cdk1"   "Smc4"   "Ckap2l" "Cks2"
</div>
```r
# Create our Seurat object and complete the initalization steps
experiment.aggregate <- CellCycleScoring(experiment.aggregate, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
```

<div class='r_output'> Warning: The following features are not present in the object: Cdca7, Cenpu,
 Dscc1, Clspn, Ung, Dtl, Exo1, E2f8, Cdc6, Cdc45, Uhrf1, not searching for symbol
 synonyms
</div>
<div class='r_output'> Warning: The following features are not present in the object: Cdk1, Ckap2l,
 Top2a, Nusap1, Cdc25c, Ndc80, Hmmr, Cenpe, Cks1brt, Ube2c, Nek2, Kif2c, Kif20b,
 Aurkb, Nuf2, Psrc1, Cdca2, Gtse1, Pimreg, Bub1, Ttk, Dlgap5, Cdc20, Ccnb2, not
 searching for symbol synonyms
</div>
```r
table(scran=experiment.aggregate$cell.cycle, seurat=experiment.aggregate$Phase)
```

<div class='r_output'>      seurat
 scran   G1  G2M    S
   G1   982 1072  699
   G2M   45   49   17
   S      9   13   10
</div>
```r
table(Idents(experiment.aggregate))
```

<div class='r_output'> 
  G2M   G1    S 
 1134 1036  726
</div>
```r
## So lets change it back to samplename
Idents(experiment.aggregate) <- "orig.ident"
table(Idents(experiment.aggregate))
```

<div class='r_output'> 
  UCD_Adj_VitE UCD_Supp_VitE  UCD_VitE_Def 
           904          1000           992
</div>
### Lets create a fake "batch" metadata (used in part 3), Here we determine UCD_Adj_VitE is from one batch and UCD_Adj_VitE/UCD_Adj_VitE are from a second battch

Here we build a new metadata variable 'batchid' which can be used to specify treatment groups.

```r
samplename = experiment.aggregate$orig.ident

batchid = rep("Batch1",length(samplename))
batchid[samplename %in% c("UCD_Adj_VitE")] = "Batch2"
names(batchid) = colnames(experiment.aggregate)

experiment.aggregate <- AddMetaData(
  object = experiment.aggregate,
  metadata = batchid,
  col.name = "batchid")

table(experiment.aggregate$batchid)
```

<div class='r_output'> 
 Batch1 Batch2 
   1992    904
</div>
### Lets spend a little time getting to know the Seurat object.

The Seurat object is the center of each single cell analysis. It stores __all__ information associated with the dataset, including data, annotations, analyses, etc. The R function slotNames can be used to view the slot names within an object.


```r
slotNames(experiment.aggregate)
```

<div class='r_output'>  [1] "assays"       "meta.data"    "active.assay" "active.ident" "graphs"      
  [6] "neighbors"    "reductions"   "project.name" "misc"         "version"     
 [11] "commands"     "tools"
</div>

```r
head(experiment.aggregate[[]])
```

<div class='r_output'>                                  orig.ident nCount_RNA nFeature_RNA
 ACTCTAATGTGGGTATG-UCD_Adj_VitE UCD_Adj_VitE       8885         2818
 AGGCTGGTCAATCACAC-UCD_Adj_VitE UCD_Adj_VitE       5019         2167
 ATGACTAGCACATGACT-UCD_Adj_VitE UCD_Adj_VitE       2208         1286
 AAGCGTCGTCTCTAAGG-UCD_Adj_VitE UCD_Adj_VitE       2795         1474
 ACATCGGGTCCATGCTC-UCD_Adj_VitE UCD_Adj_VitE       5372         2271
 ATACGGTAGTGACCAAG-UCD_Adj_VitE UCD_Adj_VitE        598          367
                                percent.mito cell.cycle      S.Score
 ACTCTAATGTGGGTATG-UCD_Adj_VitE     1.969612         G1  0.055427130
 AGGCTGGTCAATCACAC-UCD_Adj_VitE     6.216378         G1 -0.003073448
 ATGACTAGCACATGACT-UCD_Adj_VitE     6.838768         G1 -0.035830619
 AAGCGTCGTCTCTAAGG-UCD_Adj_VitE     4.221825         G1  0.004570768
 ACATCGGGTCCATGCTC-UCD_Adj_VitE     7.557707         G1  0.002626878
 ATACGGTAGTGACCAAG-UCD_Adj_VitE    11.371237         G1 -0.016286645
                                   G2M.Score Phase    old.ident batchid
 ACTCTAATGTGGGTATG-UCD_Adj_VitE  0.243198159   G2M UCD_Adj_VitE  Batch2
 AGGCTGGTCAATCACAC-UCD_Adj_VitE -0.039775669    G1 UCD_Adj_VitE  Batch2
 ATGACTAGCACATGACT-UCD_Adj_VitE  0.001092896   G2M UCD_Adj_VitE  Batch2
 AAGCGTCGTCTCTAAGG-UCD_Adj_VitE -0.041731378     S UCD_Adj_VitE  Batch2
 ACATCGGGTCCATGCTC-UCD_Adj_VitE -0.150129422     S UCD_Adj_VitE  Batch2
 ATACGGTAGTGACCAAG-UCD_Adj_VitE -0.015530630    G1 UCD_Adj_VitE  Batch2
</div>
#### Question(s)

1. What slots are empty, what slots have data?
2. What columns are available in meta.data?
3. Look up the help documentation for subset?

## Finally, save the original object, write out a tab-delimited table that could be read into excel, and view the object.

```r
write.table(as.matrix(experiment.data),"raw.datatable.txt",sep="\t",col.names=T,row.names=T)
experiment.aggregate
```

<div class='r_output'> An object of class Seurat 
 12811 features across 2896 samples within 1 assay 
 Active assay: RNA (12811 features, 0 variable features)
</div>
```r
## Original dataset in Seurat class, with no filtering
save(experiment.aggregate,file="original_seurat_object.RData")
```

## Get the next Rmd file

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/master/scrnaseq_analysis/scRNA_Workshop-PART2.Rmd", "scRNA_Workshop-PART2.Rmd")
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
 [1] biomaRt_2.44.0 Seurat_3.1.5  
 
 loaded via a namespace (and not attached):
   [1] ggbeeswarm_0.6.0            Rtsne_0.15                 
   [3] colorspace_1.4-1            ellipsis_0.3.1             
   [5] ggridges_0.5.2              XVector_0.28.0             
   [7] BiocNeighbors_1.6.0         GenomicRanges_1.40.0       
   [9] leiden_0.3.3                listenv_0.8.0              
  [11] ggrepel_0.8.2               bit64_0.9-7                
  [13] AnnotationDbi_1.50.0        codetools_0.2-16           
  [15] splines_4.0.0               scater_1.16.0              
  [17] knitr_1.28                  jsonlite_1.6.1             
  [19] ica_1.0-2                   cluster_2.1.0              
  [21] dbplyr_1.4.3                png_0.1-7                  
  [23] uwot_0.1.8                  sctransform_0.2.1          
  [25] BiocManager_1.30.10         compiler_4.0.0             
  [27] httr_1.4.1                  dqrng_0.2.1                
  [29] assertthat_0.2.1            Matrix_1.2-18              
  [31] lazyeval_0.2.2              limma_3.44.1               
  [33] BiocSingular_1.4.0          htmltools_0.4.0            
  [35] prettyunits_1.1.1           tools_4.0.0                
  [37] rsvd_1.0.3                  igraph_1.2.5               
  [39] gtable_0.3.0                glue_1.4.1                 
  [41] GenomeInfoDbData_1.2.3      RANN_2.6.1                 
  [43] reshape2_1.4.4              dplyr_0.8.5                
  [45] rappdirs_0.3.1              Rcpp_1.0.4.6               
  [47] Biobase_2.48.0              vctrs_0.3.0                
  [49] ape_5.3                     nlme_3.1-147               
  [51] DelayedMatrixStats_1.10.0   lmtest_0.9-37              
  [53] xfun_0.13                   stringr_1.4.0              
  [55] globals_0.12.5              lifecycle_0.2.0            
  [57] irlba_2.3.3                 renv_0.10.0                
  [59] statmod_1.4.34              XML_3.99-0.3               
  [61] future_1.17.0               edgeR_3.30.0               
  [63] zlibbioc_1.34.0             MASS_7.3-51.5              
  [65] zoo_1.8-8                   scales_1.1.1               
  [67] hms_0.5.3                   parallel_4.0.0             
  [69] SummarizedExperiment_1.18.1 RColorBrewer_1.1-2         
  [71] SingleCellExperiment_1.10.1 yaml_2.2.1                 
  [73] curl_4.3                    memoise_1.1.0              
  [75] reticulate_1.15             pbapply_1.4-2              
  [77] gridExtra_2.3               ggplot2_3.3.0              
  [79] stringi_1.4.6               RSQLite_2.2.0              
  [81] S4Vectors_0.26.1            scran_1.16.0               
  [83] BiocGenerics_0.34.0         BiocParallel_1.22.0        
  [85] GenomeInfoDb_1.24.0         matrixStats_0.56.0         
  [87] rlang_0.4.6                 pkgconfig_2.0.3            
  [89] bitops_1.0-6                evaluate_0.14              
  [91] lattice_0.20-41             ROCR_1.0-11                
  [93] purrr_0.3.4                 patchwork_1.0.0            
  [95] htmlwidgets_1.5.1           cowplot_1.0.0              
  [97] bit_1.1-15.2                tidyselect_1.1.0           
  [99] RcppAnnoy_0.0.16            plyr_1.8.6                 
 [101] magrittr_1.5                R6_2.4.1                   
 [103] IRanges_2.22.1              DelayedArray_0.14.0        
 [105] DBI_1.1.0                   pillar_1.4.4               
 [107] fitdistrplus_1.1-1          survival_3.1-12            
 [109] RCurl_1.98-1.2              tibble_3.0.1               
 [111] future.apply_1.5.0          tsne_0.1-3                 
 [113] crayon_1.3.4                KernSmooth_2.23-16         
 [115] BiocFileCache_1.12.0        plotly_4.9.2.1             
 [117] rmarkdown_2.1               viridis_0.5.1              
 [119] progress_1.2.2              locfit_1.5-9.4             
 [121] grid_4.0.0                  data.table_1.12.8          
 [123] blob_1.2.1                  digest_0.6.25              
 [125] tidyr_1.0.3                 openssl_1.4.1              
 [127] stats4_4.0.0                munsell_0.5.0              
 [129] beeswarm_0.2.3              viridisLite_0.3.0          
 [131] vipor_0.4.5                 askpass_1.1
</div>