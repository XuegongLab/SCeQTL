#' Normalize gene expression matrix
#'
#' Normalize the gene expression matrix using DESeq method. You could find detailed description in paper "Differential expression analysis for sequence count data"
#'
#' @param gene A gene expression matrix where each row is one gene and each column is one sample
#' @return 
#' A gene expression matrix after nomalization. Genes that not express in every samples will be filtered.
#' @examples
#' gene <- normalize(gene)
normalize <- function(gene) {
  geo_mean <- function(data) {
    log_data <- log(data)
    gm <- exp(mean(log_data[is.finite(log_data)]))
    return(gm)
  }
  data = gene
  rowsum = apply(data, 1, sum)
  data = data[rowsum!=0,]

  gene.size = apply(data, 1, geo_mean)
  data1 = data/gene.size
  librarysize = apply(data1, 2, function(x) median(x[x!=0]))

  data = floor(data/librarysize)

  rowsum = apply(data, 1, sum)
  data_filter_0_gene = data[rowsum!=0,]

  return(data_filter_0_gene)
}

