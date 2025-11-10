---
title: "GREAT"
author: Nathan Harmston
date: '05 June 2024'
format:
  html: 
    keep-md: true
    embed-resources: true
    df-print: kable
    toc: true
    toc-depth: 3
    code-fold: true
    number-sections: true
    smooth-scroll: true
    code-tools: true
    code-line-numbers: true
  gfm: 
    df-print: kable
    toc: true 
    toc-depth: 3
    number-sections: true
---



# GREAT analysis of deaf1 peaks 



::: {.cell}

```{.r .cell-code}
library(rtracklayer)
```

::: {.cell-output .cell-output-stderr}
```
Loading required package: GenomicRanges
```
:::

::: {.cell-output .cell-output-stderr}
```
Loading required package: stats4
```
:::

::: {.cell-output .cell-output-stderr}
```
Loading required package: BiocGenerics
```
:::

::: {.cell-output .cell-output-stderr}
```

Attaching package: 'BiocGenerics'
```
:::

::: {.cell-output .cell-output-stderr}
```
The following objects are masked from 'package:stats':

    IQR, mad, sd, var, xtabs
```
:::

::: {.cell-output .cell-output-stderr}
```
The following objects are masked from 'package:base':

    anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    table, tapply, union, unique, unsplit, which.max, which.min
```
:::

::: {.cell-output .cell-output-stderr}
```
Loading required package: S4Vectors
```
:::

::: {.cell-output .cell-output-stderr}
```
Warning: package 'S4Vectors' was built under R version 4.3.2
```
:::

::: {.cell-output .cell-output-stderr}
```

Attaching package: 'S4Vectors'
```
:::

::: {.cell-output .cell-output-stderr}
```
The following object is masked from 'package:utils':

    findMatches
```
:::

::: {.cell-output .cell-output-stderr}
```
The following objects are masked from 'package:base':

    expand.grid, I, unname
```
:::

::: {.cell-output .cell-output-stderr}
```
Loading required package: IRanges
```
:::

::: {.cell-output .cell-output-stderr}
```
Loading required package: GenomeInfoDb
```
:::

::: {.cell-output .cell-output-stderr}
```
Warning: package 'GenomeInfoDb' was built under R version 4.3.3
```
:::

```{.r .cell-code}
library(BSgenome.Mmusculus.UCSC.mm10)
```

::: {.cell-output .cell-output-stderr}
```
Loading required package: BSgenome
```
:::

::: {.cell-output .cell-output-stderr}
```
Warning: package 'BSgenome' was built under R version 4.3.2
```
:::

::: {.cell-output .cell-output-stderr}
```
Loading required package: Biostrings
```
:::

::: {.cell-output .cell-output-stderr}
```
Warning: package 'Biostrings' was built under R version 4.3.3
```
:::

::: {.cell-output .cell-output-stderr}
```
Loading required package: XVector
```
:::

::: {.cell-output .cell-output-stderr}
```

Attaching package: 'Biostrings'
```
:::

::: {.cell-output .cell-output-stderr}
```
The following object is masked from 'package:base':

    strsplit
```
:::

::: {.cell-output .cell-output-stderr}
```
Loading required package: BiocIO
```
:::

::: {.cell-output .cell-output-stderr}
```

Attaching package: 'BiocIO'
```
:::

::: {.cell-output .cell-output-stderr}
```
The following object is masked from 'package:rtracklayer':

    FileForFormat
```
:::

```{.r .cell-code}
deaf1.idr = read.delim("../peaks/deaf1-idr", header=FALSE)
#deaf1.filtered.idr = deaf1.idr[deaf1.idr[,5]>=540,] 
#deaf1.filtered.idr = deaf1.filtered.idr[!deaf1.filtered.idr[,1] %in% c("chr4_GL456216_random", "chr4_JH584295_random"), ] 
 
deaf1.gr = GRanges(deaf1.idr[,1], IRanges(deaf1.idr$V2+1, deaf1.idr$V3), idr=deaf1.idr[,5], summit=deaf1.idr[,10])

rep1 = import("../peaks/deaf1_rep1_peaks.narrowPeak")
rep2 = import("../peaks/deaf1_rep2_peaks.narrowPeak")

idr.ol1 = findOverlaps(deaf1.gr, rep1)
idr.ol2 = findOverlaps(deaf1.gr, rep2)

deaf1.gr$rep1_pvalue[queryHits(idr.ol1)] =  rep1$pValue[subjectHits(idr.ol1)]
deaf1.gr$rep2_pvalue[queryHits(idr.ol2)] =  rep2$pValue[subjectHits(idr.ol2)]

deaf1.gr$rep1_qvalue[queryHits(idr.ol1)] =  rep1$qValue[subjectHits(idr.ol1)]
deaf1.gr$rep2_qvalue[queryHits(idr.ol2)] =  rep1$qValue[subjectHits(idr.ol2)]
  
sum( deaf1.gr$rep1_qvalue > -log10(0.05) & deaf1.gr$rep2_qvalue > -log10(0.05))
```

::: {.cell-output .cell-output-stdout}
```
[1] 22884
```
:::

```{.r .cell-code}
deaf1.filtered.gr = deaf1.gr[deaf1.gr$idr >= 540]

#sum( deaf1.filtered.gr$rep1_qvalue > -log10(0.01) & deaf1.filtered.gr$rep2_qvalue > -log10(0.01))

deaf1.filtered.gr = deaf1.filtered.gr[!seqnames(deaf1.filtered.gr) %in% c("chr4_GL456216_random", "chr4_JH584295_random"), ] 

names(deaf1.filtered.gr) = paste("peak", 1:length(deaf1.filtered.gr), sep=":")
names(deaf1.filtered.gr) = paste("peak", 1:length(deaf1.filtered.gr), sep=":")

export.bed(deaf1.filtered.gr, "../peaks/deaf1_idr.bed")
write.csv(as.data.frame(deaf1.filtered.gr), "../peaks/deaf1_idr.csv", quote=FALSE)


deaf1.filtered.gr$new_start = start(deaf1.filtered.gr) + deaf1.filtered.gr$summit


deaf1.gr = GRanges(seqnames(deaf1.filtered.gr), IRanges(deaf1.filtered.gr$new_start-250, deaf1.filtered.gr$new_start+250))

peak.sequences = getSeq(Mmusculus, deaf1.gr)
names(peak.sequences) = paste("peak", 1:length(peak.sequences), sep=":")
writeXStringSet(peak.sequences, "../peaks/deaf1_idr.fasta")
```
:::

::: {.cell}

```{.bash .cell-code}

```
:::

::: {.cell}

```{.r .cell-code}
sessionInfo()
```

::: {.cell-output .cell-output-stdout}
```
R version 4.3.1 (2023-06-16)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/London
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] BSgenome.Mmusculus.UCSC.mm10_1.4.3 BSgenome_1.70.2                   
 [3] BiocIO_1.12.0                      Biostrings_2.70.3                 
 [5] XVector_0.42.0                     rtracklayer_1.62.0                
 [7] GenomicRanges_1.54.1               GenomeInfoDb_1.38.8               
 [9] IRanges_2.36.0                     S4Vectors_0.40.2                  
[11] BiocGenerics_0.48.1               

loaded via a namespace (and not attached):
 [1] Matrix_1.6-5                jsonlite_1.8.8             
 [3] compiler_4.3.1              rjson_0.2.21               
 [5] crayon_1.5.2                SummarizedExperiment_1.32.0
 [7] Biobase_2.62.0              Rsamtools_2.18.0           
 [9] bitops_1.0-7                GenomicAlignments_1.38.2   
[11] parallel_4.3.1              BiocParallel_1.36.0        
[13] yaml_2.3.8                  fastmap_1.2.0              
[15] lattice_0.22-6              S4Arrays_1.2.1             
[17] knitr_1.47                  htmlwidgets_1.6.4          
[19] XML_3.99-0.16.1             DelayedArray_0.28.0        
[21] MatrixGenerics_1.14.0       GenomeInfoDbData_1.2.11    
[23] rlang_1.1.3                 xfun_0.44                  
[25] SparseArray_1.2.4           cli_3.6.2                  
[27] zlibbioc_1.48.2             grid_4.3.1                 
[29] digest_0.6.35               rstudioapi_0.16.0          
[31] evaluate_0.23               codetools_0.2-20           
[33] abind_1.4-5                 RCurl_1.98-1.14            
[35] restfulr_0.0.15             rmarkdown_2.27             
[37] matrixStats_1.3.0           tools_4.3.1                
[39] htmltools_0.5.8.1          
```
:::
:::
