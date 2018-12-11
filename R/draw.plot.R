drawhistplot <- function(sample.data,removeZero=TRUE){
  sample.data.no.zero = sample.data
  if(removeZero){
    sample.data.no.zero = sample.data[sample.data$expression!=0,]
  }
  ggplot(sample.data.no.zero, aes(expression, fill = snp))+geom_histogram()+facet_grid(snp ~ ., margins=TRUE, scales="free_y")
}

drawboxplot <- function(sample.data,removeZero=TRUE){
  sample.data.no.zero = sample.data
  if(removeZero){
    sample.data.no.zero = sample.data[sample.data$expression!=0,]
  }
  ggplot(sample.data.no.zero, aes(factor(snp),expression))+geom_boxplot()
}

drawviolinplot <- function(sample.data,removeZero=TRUE){
  sample.data.no.zero = sample.data
  if(removeZero){
    sample.data.no.zero = sample.data[sample.data$expression!=0,]
  }
  ggplot(sample.data.no.zero, aes(factor(snp),expression))+geom_violin()
}

checkdist <- function(gene, n=10){
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

