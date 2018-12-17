# SCeQTL <img src="https://github.com/sunfenghao2017/logo/blob/master/logo.png" align="right" height =  50 width= 150/>

[![Build Status](https://travis-ci.org/r-lib/devtools.svg?branch=master)](https://travis-ci.org/r-lib/SCeQTL)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/r-lib/SCeQTL?branch=master&svg=true)](https://ci.appveyor.com/project/hadley/SCeQTL)
[![Coverage Status](https://codecov.io/github/r-lib/SCeQTL/coverage.svg?branch=master)](https://codecov.io/github/r-lib/SCeQTL?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/SCeQTL)](https://cran.r-project.org/package=SCeQTL)

SCeQTL is an R package that uses zero-inflated negative binomial regression to do eQTL analysis on single-cell data. It can distinguish two type of gene-expression differences among different genotype groups. It’s more suitable to use SCeQTL to identify eQTLs from single-cell data. It can also be used for finding gene expression variations associated with other grouping factors like cell lineages. Following is the detail information and usage of this program. You can also found them in READ.ME. [R
Packages](http://r-pkgs.had.co.nz/) is a book based around this workflow.
## Authors
Hu Yue <y-hu10@qq.com>

## Version

For **Windows**:[SCeQTL.zip](http://)

For **Linux and Mac**:[SCeQTL.tar.gz](http://)

## Installation

```r
#You can download the R package from here and run the command in R:

install.packages(<source_code_file>, repos=NULL, type="source")


```
## Requirement

### Prerequisite(R-packages): 
       qvalue
       ggplot2
       doParallel
       pscl

## Function:

The function of this program is to test the correlation between gene splicing expression and gene variants.

**input**: Rdata which include gene expression matrix and SNV matrix

**output**: P value or q value of gene-SNV pairs

## Usage
```r
SCeQTL:: cal.pvalue <- function(gene, SNV, thread = 8, remove_outlier = TRUE,EM = TRUE, dist = 'negbin', type = 0)

#Calculatre and return P value of gene-SNV pairs.
```

* `thread`: number of threads that are used.

* `remove_outlier`: whether program should remove samples whose expression level are far away from the others(>median+4*mad)

* `EM`: use EM or BGFS to optimaize

* `Dist`: distribution assumption, could be negative binomial or Poisson

* `type`: which kind of difference you are interested in. type 0 means non-zero part difference, type 1 means zero ratio difference, type 2 means at least one or non-zero part or zero ratio difference

 
```r
SCeQTL::checkdist(gene, n=10)

#The function samples n gene expressions, draw normalized QQ-plot to compare real gene distribution with fitted gene distribution. The function is used for checking whether non-zero part of the data fit negative binomial distribution well.
```
 
```r
SCeQTL::check.sample <- function(sample_gene, sample_SNV, plottype='boxplot', removeZero = TRUE)

#Print statistic and visualize one gene-SNV pair.
```
### Input format:


**Rdata** which stored the **gene expression matrix** and **SNV matrix**.
    
---
**1.** * `Gene expression matrix`   should be named as “gene”, where one row indicate one gene and one column indicate one sample.
 
**2.** * `SNV matrix`                 should be named as “snp”, where one row indicate one variant and one column indicate one sample.
### Test data:

You can find test data here.

## Code of conduct

Please note that the SCeQTL project is released with a [Contributor Code of Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.
