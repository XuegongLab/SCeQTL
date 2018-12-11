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

checkzeroratio <- function(sample.data){
  cat("zero ratio of group 0",sum(sample.data[sample.data$snp==0,1]==0)/sum(sample.data[sample.data$snp==0,1]>=0),"\n")
  cat("zero ratio of group 1",sum(sample.data[sample.data$snp==1,1]==0)/sum(sample.data[sample.data$snp==1,1]>=0),"\n")
  cat("zero ratio of group 2",sum(sample.data[sample.data$snp==2,1]==0)/sum(sample.data[sample.data$snp==2,1]>=0),"\n")
}

check.sample <- function(sample_gene, sample_snp, plottype='boxplot', removeZero = TRUE){
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
