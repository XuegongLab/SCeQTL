#' Print statistic and visualize one gene-SNV pair.
#'
#' @param sample_gene gene expression vector that each element in it is the gene expression level in one sample
#' @param sample_snp snp vector that each element in it is the genotype in one sample
#' @param plottype an optional parameter used to choose the type of plot, candidates are \code{boxplot}, \code{violin}, \code{histplot}. \code{boxplot} is used if not specified.
#' @param removeZero an optional parameter used to decide if you'd like to exclude non-expressed samples when plot the sample. Setting to \code{TRUE} helps to check non-zero part difference while setting to \code{FALSE} helps to check zero ratio difference.\code{TRUE} is used if not specified.
#' @return 
#' NONE.
#' output the mean of each group, non-zero part mean of each group, non-zero part standard error of each group to the terminal and a picture show the relationship between gene and snp just like picture below. 
#' @note 
#' Box Plot \figure{boxplot.png}
#' @examples
#' check.sample(gene[1,], snp[1,])
#' check.sample(gene[1,], snp[1,], 'histplot')
#' check.sample(gene[1,], snp[1,], 'histplot', FALSE)
check.sample <- function(sample_gene, sample_snp, plottype='boxplot', removeZero = TRUE){
  checkzeroratio <- function(sample.data){
    cat("zero ratio of group 0",sum(sample.data[sample.data$snp==0,1]==0)/sum(sample.data[sample.data$snp==0,1]>=0),"\n")
    cat("zero ratio of group 1",sum(sample.data[sample.data$snp==1,1]==0)/sum(sample.data[sample.data$snp==1,1]>=0),"\n")
    cat("zero ratio of group 2",sum(sample.data[sample.data$snp==2,1]==0)/sum(sample.data[sample.data$snp==2,1]>=0),"\n")
  }
  
  checkmeanvar <- function(sample.data){
    cat("mean of group 0: ",mean(sample.data[sample.data$snp==0,1]),"\n")
    cat("mean of group 1: ",mean(sample.data[sample.data$snp==1,1]),"\n")
    cat("mean of group 2: ",mean(sample.data[sample.data$snp==2,1]),"\n")
    sample.data.no.zero = sample.data[sample.data$expression!=0,]
    cat("non-zero part mean of group 0: ",mean(sample.data.no.zero[sample.data.no.zero$snp==0,1]),"\n")
    cat("non-zero part mean of group 1: ",mean(sample.data.no.zero[sample.data.no.zero$snp==1,1]),"\n")
    cat("non-zero part mean of group 2: ",mean(sample.data.no.zero[sample.data.no.zero$snp==2,1]),"\n")
    cat("non-zero part standard error of group 0: ",sd(sample.data.no.zero[sample.data.no.zero$snp==0,1]),"\n")
    cat("non-zero part standard error of group 1: ",sd(sample.data.no.zero[sample.data.no.zero$snp==1,1]),"\n")
    cat("non-zero part standard error of group 2: ",sd(sample.data.no.zero[sample.data.no.zero$snp==2,1]),"\n")
  }
  
  drawhistplot <- function(sample.data,removeZero=TRUE){
    sample.data.no.zero = sample.data
    if(removeZero){
      sample.data.no.zero = sample.data[sample.data$expression!=0,]
    }
    snp = NULL
    ggplot(sample.data.no.zero, aes(expression, fill = snp))+geom_histogram()+facet_grid(snp ~ ., margins=TRUE, scales="free_y")
  }
  
  drawboxplot <- function(sample.data,removeZero=TRUE){
    sample.data.no.zero = sample.data
    if(removeZero){
      sample.data.no.zero = sample.data[sample.data$expression!=0,]
    }
    snp = NULL
    ggplot(sample.data.no.zero, aes(factor(snp),expression))+geom_boxplot()
  }
  
  drawviolinplot <- function(sample.data,removeZero=TRUE){
    sample.data.no.zero = sample.data
    if(removeZero){
      sample.data.no.zero = sample.data[sample.data$expression!=0,]
    }
    snp = NULL
    ggplot(sample.data.no.zero, aes(factor(snp),expression))+geom_violin()
  }
  
  removeoutlier <- function(sample.data){
    sample.data.no.zero = sample.data[sample.data$expression!=0,1]
    med = median(sample.data.no.zero)
    Mad = mad(sample.data.no.zero)
    return(sample.data[sample.data$expression<med+4*Mad,])
  }
  
  sample.data = data.frame(unlist(sample_gene), unlist(sample_snp))
  colnames(sample.data) = c('expression','snp')
  sample.data = removeoutlier(sample.data)

  checkzeroratio(sample.data)
  checkmeanvar(sample.data)
  
  if(plottype=='violin')
    drawviolinplot(sample.data,removeZero = removeZero)
  else if(plottype=='boxplot')
    drawboxplot(sample.data, removeZero = removeZero)
  else if(plottype=='histplot')
    drawhistplot(sample.data, removeZero = removeZero)
}
