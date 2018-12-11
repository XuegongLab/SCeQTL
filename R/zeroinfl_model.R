zeroinfl_model <- function(sample_gene, sample_snp, remove_outlier=TRUE, dist='negbin', EM = TRUE, type = 0){
  sample.data = data.frame(unlist(sample_gene), unlist(sample_snp))
  colnames(sample.data) = c('expression','snp')
  if(remove_outlier)
    sample.data = removeoutlier(sample.data)
  m1 <- try(zeroinfl(expression ~ snp, data = sample.data, dist = dist, EM = EM), silent=TRUE)
  if(class(m1)=="try-error")
    return(NA)

  if(type==0){
    m0 <- zeroinfl(expression ~ snp|1, data = sample.data, dist = dist, EM = EM)
    .df <- 1
  }
  else if(type==1){
    m0 <- zeroinfl(expression ~ 1|snp, data = sample.data, dist = dist, EM = EM)
    .df <- 1
  }
  else{
    m0 <- zeroinfl(expression ~ 1, data = sample.data, dist = dist, EM = EM)
    .df <- 2
  }
  return(pchisq(2 * (logLik(m1) - logLik(m0)), df = .df, lower.tail=FALSE))
}
