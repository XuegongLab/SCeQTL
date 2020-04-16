#' Calculatre and return P value of gene-SNV pairs.
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
#' }
#' @seealso 
#' Before calculating p-value, you may want to check whether Poisson or negative bionomial is better using \link{checkdist}. And you could normalize the gene expression matrix by \link{normalize}. You could also calculate q value by \link{cal.qvalue}. After you get significant result, you could check it by \link{check.sample}.
#' @examples
#' pvalue <- cal.pvalue(gene, snp)
#' pvalue <- cal.pvalue(gene, snp, thread = 16)
#' pvalue <- cal.pvalue(gene, snp, remove_outlier = FALSE)
#' pvalue <- cal.pvalue(gene, snp, EM = FALSE)
#' pvalue <- cal.pvalue(gene, snp, dist = 'poisson')
#' pvalue <- cal.pvalue(gene, snp, type = 2)
cal.pvalue <- function(gene, snp, thread = 8, remove_outlier = TRUE,EM = TRUE, dist = 'negbin', type = 0){
    removeoutlier <- function(sample.data){
        sample.data.no.zero = sample.data[sample.data$expression!=0,1]
        med = median(sample.data.no.zero)
        Mad = mad(sample.data.no.zero)
        return(sample.data[sample.data$expression<med+4*Mad,])
    }

    zeroinfl_model <- function(sample_gene, sample_snp, remove_outlier=TRUE, dist='negbin', EM = TRUE, type = 0){
        sample.data = data.frame(unlist(sample_gene), unlist(sample_snp))
        colnames(sample.data) = c('expression','snp')
        if(remove_outlier)
          sample.data = removeoutlier(sample.data)
        m1 <- try(zeroinfl(expression ~ snp, data = sample.data, dist = dist, EM = EM), silent=TRUE)
        if(class(m1)=="try-error")
          return(NA)

        if(type==0){
          m0 <- try(zeroinfl(expression ~ snp|1, data = sample.data, dist = dist, EM = EM))
          if(class(pvalue)=="try-error"){
            message("\twarning: no valid set of coefficients has been found\n")
            return(NA)
          }
          .df <- 1
        }
        else if(type==1){
          m0 <- try(zeroinfl(expression ~ 1|snp, data = sample.data, dist = dist, EM = EM))
          if(class(pvalue)=="try-error"){
            message("\twarning: no valid set of coefficients has been found\n")
            return(NA)
          }
          .df <- 1
        }
        else{
          m0 <- try(zeroinfl(expression ~ 1, data = sample.data, dist = dist, EM = EM))
          if(class(pvalue)=="try-error"){
            message("\twarning: no valid set of coefficients has been found\n")
            return(NA)
          }
          .df <- 2
        }
        pvalue = try(pchisq(2 * (logLik(m1) - logLik(m0)), df = .df, lower.tail=FALSE))
        if(class(pvalue)=="try-error"){
          message("\twarning: no valid set of coefficients has been found\n")
          return(NA)
        }
        else
          return(pvalue)
    }
    if (type == 0) 
        message("Identyfing non-zero part difference...\n")
    else if (type == 1) 
        message("Identyfing zero ratio difference...\n")
    else 
        message("Identyfing non-zero part or/and zero ratio difference...\n")
    registerDoParallel(thread)
    countzero <- rowSums(gene != 0)
    gene <- gene[countzero > 3,]
    removed <- labels(countzero)[countzero <= 3]
    if (length(removed) != 0)
        for(i in 1:length(removed))
            message(paste0("Gene: ", removed[i], " was removed because of excess zero.\n"))
    else 
        message('No gene was removed because of excess zero.\n')
    if (nrow(gene) == 0) stop('No gene to be tested.')
    gene.count = dim(gene)[1]
    snp.count = dim(snp)[1]
    pvalue = list()
    j = 0
    message("Start calculating p value...\n")
    for(i in 1:gene.count){
        message(paste0("calculate pvalue for gene: ", i, "\n"))
        result = foreach(j=1:snp.count) %dopar% {zeroinfl_model(gene[i,],snp[j,],remove_outlier=remove_outlier,dist=dist,EM=EM,type=type)}
        pvalue = rbind(pvalue,result)
    }
    gene.name = rep(row.names(gene), each = snp.count)
    if(is.null(row.names(snp))){
        snp.raw.name = 1:snp.count
    }else
        snp.raw.name = row.names(snp)
    snp.name = list()
    for(i in 1:gene.count)
        snp.name = c(snp.name, snp.raw.name)
    result = data.frame(gene.name, unlist(snp.name), unlist(pvalue))
    colnames(result) <- c("gene","snp","pvalue")
    return(result)
}
