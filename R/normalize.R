geo_mean <- function(data) {
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}

normalize <- function(data) {
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

