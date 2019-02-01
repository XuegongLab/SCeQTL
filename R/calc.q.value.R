#' Calculatre and return Q value of gene-SNV pairs.
#'
#' 
#'
#' @param gene A gene expression matrix where each row is one gene and each column is one sample
#' @param snp A genotype matrix where each row is one variant and each column is one sample
#' @param thread Number of threads that are used. DEFAULT option is 8.
#' @param remove_outlier whether program should remove samples whose expression level are far away from the others(>median+4*mad). DEFAULT option is TRUE.
#' @param EM use EM(TRUE) or BGFS(FALSE) to optimaize. Default option is TRUE
#' @param dist distribution assumption, could be negative binomial(negbin) or Poisson(poisson). Default option is negative bionomial.
#' @param type which kind of difference you are interested in. type 0 means non-zero part difference, type 1 means zero ratio difference, type 2 means at least one or non-zero part or zero ratio difference. Default option is 0.
#' @return 
#' A data frame containing the SCeQTL result. Rows are gene-snp pair and columns contain the following
#' \itemize{
#'  \item gene the gene name
#'  \item snp the snp name
#'  \item pvalue the p-value of the gene-snp pair
#'  \item qvalue the q-value of the gene-snp pair
#' }
#' @seealso 
#' Before calculating q-value, you may want to check whether Poisson or negative bionomial is better using \link{checkdist}. And you could normalize the gene expression matrix by \link{normalize}. You could also calculate p value by \link{cal.pvalue}. After you get significant result, you could check it by \link{check.sample}.
#' @examples
#' qvalue <- cal.qvalue(gene, snp)
#' qvalue <- cal.qvalue(gene, snp, thread = 16)
#' qvalue <- cal.qvalue(gene, snp, remove_outlier = FALSE)
#' qvalue <- cal.qvalue(gene, snp, EM = FALSE)
#' qvalue <- cal.qvalue(gene, snp, dist = 'poisson')
#' qvalue <- cal.qvalue(gene, snp, type = 2)
cal.qvalue <- function(gene, snp, thread = 8, remove_outlier = TRUE,EM = TRUE, dist = 'negbin', type = 0){
  result = cal.pvalue(gene,snp,thread,remove_outlier,EM,dist,type)
  result$qvalue = qvalue(result$pvalue)$qvalues
  return(result)
}