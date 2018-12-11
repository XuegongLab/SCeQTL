find.eqtl <- function(){
  load('gene.snp.Rdata')
  gene <- normalize(gene)
  checkdist(gene)
  result <- cal.pvalue(gene,snp)
  return(result)
}

cal.pvalue <- function(gene, snp, thread = 8, remove_outlier = TRUE,EM = TRUE, dist = 'negbin', type = 0){
  registerDoParallel(thread)
  gene.count = dim(gene)[1]
  snp.count = dim(snp)[1]
  pvalue = list()
  message("starting calculating p value...\n")
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
    snp.name = rbind(snp.name, snp.raw.name)
  result = data.frame(gene.name, unlist(snp.name), unlist(pvalue))
  colnames(result) <- c("gene","snp","pvalue")
  return(result)
}

cal.qvalue <- function(gene, snp, thread = 8, remove_outlier = TRUE,EM = TRUE, dist = 'negbin', type = 0){
  result = cal.pvalue(gene,snp,thread,remove_outlier,EM,dist,type)
  result$qvalue = qvalue(result$pvalue)
  return(result)
}
