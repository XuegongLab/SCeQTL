removeoutlier <- function(sample.data){
  sample.data.no.zero = sample.data[sample.data$expression!=0,1]
  med = median(sample.data.no.zero)
  Mad = mad(sample.data.no.zero)
  return(sample.data[sample.data$expression<med+4*Mad,])
}
