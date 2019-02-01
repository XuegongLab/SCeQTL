#' Check gene non-zero part distribution
#'
#' The function samples n gene expressions, draw normalized QQ-plot to compare real gene distribution with fitted gene distribution. The function is used for checking whether non-zero part of the data fit negative binomial distribution well.
#'
#' @param gene The gene expression matrix you'd like to check
#' @param n The num of genes that will be randomly picked, default value is 10
#' @return 
#' NONE, a QQ-plot will be shown
#' @note 
#' This function may fail since It's possible that the random picked gene can't be fit to negative bionomial distribution, all zero value for example.
#' QQ-plot \figure{qqplot.png}
#' @examples
#' checkdist(gene)
#' checkdist(gene, 5)
checkdist <- function(gene, n=10){
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
  gene.count <- dim(gene)[1]
  sample.index <- sample(1:gene.count, n)
  fx <- c()
  fy <- c()
  index <- c()
  for(i in sample.index){
    sample.gene = gene[i,]
    sample.data = data.frame(unlist(sample.gene), 1:dim(sample.gene)[2])
    colnames(sample.data) <- c("expression","index")
    sample.data = removeoutlier(sample.data)
    m0 <- zeroinfl(expression~1, data=sample.data, dist="negbin")
    x <- unlist(sample.data$expression)
    y <- unlist(rnbinom(n=length(sample.data$expression),size=m0$theta, mu=exp(m0$coefficients$count)))
    maxx <- max(x)
    x <- x/maxx*1.25
    y <- y/maxx*1.25
    x <- x[x!=0]
    y <- y[y!=0]
    sx <- sort(x); sy <- sort(y)
    lenx <- length(sx)
    leny <- length(sy)
    if (leny < lenx)sx <- approx(1L:lenx, sx, n = leny)$y
    if (leny > lenx)sy <- approx(1L:leny, sy, n = lenx)$y
    fx <- c(fx,sx)
    fy <- c(fy,sy)
    index <- c(index, unlist(replicate(n = length(x),expr = i)))
  }
  qqplot.data <- data.frame(fx,fy,index)
  colnames(qqplot.data) <- c("real distribution", "fitted distribution", "index")
  ggplot(data = qqplot.data,aes(x = fx,y=fy),grouping=index)+xlim(0,1)+ylim(0,1)+geom_path()
}