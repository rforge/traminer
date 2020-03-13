## function for comparing sets of sequences by computing LRT and tail probability
## Tim Liao, University of Illionis, and Anette Fasang, Humboldt University
## version 1.0, March 2020

seqBIC <- function(seqdata, seqdata2=NULL, group=NULL, subgroup=NULL,
  s=100, seed=36963, squared="LRTonly", method, ...)
{
  return(seqLRT(seqdata, seqdata2, group, subgroup, s=s, seed=seed,
         stat="BIC", squared=squared, method, ...))
}


seqLRT <- function(seqdata, seqdata2=NULL, group=NULL, subgroup=NULL,
  s=100, seed=36963, stat="LRT", squared="LRTonly", method, ...)
{
  #require("gtools")
  #require("TraMineR")


  if (is.null(seqdata2) & is.null(group)){
    stop("'seqdata2' and 'group' cannot both be NULL!")
  }
  if (!is.null(subgroup) & is.null(group)){
    stop("'subgroup' not NULL while 'group' is NULL!")
  }

  if (is.logical(squared))
    LRTpow <- ifelse(squared, 2, 1)
  else {
    if (squared != "LRTonly") stop("squared must be logical or 'LRTonly'")
    LRTpow <- 2
    squared <- FALSE
  }

  is1.stslist <- inherits(seqdata,"stslist")
  is2.stslist <- inherits(seqdata2,"stslist")

  if (is.list(seqdata) & !is1.stslist) {
    if (is2.stslist | length(seqdata) != length(seqdata2))
      stop("When 'seqdata' is a list, seqdata2 must be a list of same length")
    else {
      l <- length(seqdata)
      test <- FALSE
      i <- 1
      while (i <= l) {
        if(!inherits(seqdata[[i]], "stslist") | !inherits(seqdata2[[i]], "stslist"))
           stop("At least one element of the seqdata lists is not a state sequence object!")
        i=i+1
      }
    }
  }
  else if (!is1.stslist) {
    stop("If not a list, 'seqdata' must be a state sequence object (stslist) created with seqdef()")
  }
  else if (!is.null(seqdata2) & !is2.stslist) {
    stop("If not a list, 'seqdata2' must be a state sequence object (stslist) created with seqdef()")
  }

  if (!stat %in% c("LRT","BIC","all"))
    stop("Bad stat value, must be one of 'LRT', 'BIC', or 'all'")

  if (stat=="all") {
    is.LRT <- is.BIC <- TRUE
  }
  else{
    is.LRT <- "LRT" %in% stat
    is.BIC <- "BIC" %in% stat
  }

  if (!is1.stslist){
    seq1 <- seqdata
    seq2 <- seqdata2
  }
  else if (is1.stslist & !is.null(seqdata2)) {
    seq1 <- list(seqdata)
    seq2 <- list(seqdata2)
  }
  else if (is1.stslist & is.null(seqdata2)) {
    ## suppress cases with NA group values
    gvar <- as.vector(group)
    if (!is.null(subgroup)){
      sgvar <- as.vector(subgroup)
      inotna <- which(!is.na(gvar) & !is.na(sgvar))
      sgvar <- sgvar[inotna]
      sgvar <- factor(sgvar)
      lev.sg <- levels(sgvar)
    }
    else {
      inotna <- which(!is.na(gvar))
    }
    ########
    ina <- nrow(seqdata) - length(inotna)
    if(ina > 0)
      message("[!!] ", ina, " sequences removed because of NA values of the grouping variable(s)\n")
    ##########
    gvar <- gvar[inotna]
    gvar <- factor(gvar)
    lev.g <- levels(gvar)
    if (length(lev.g) > 2)
      stop("Currently seqLRT supports only 2 groups!")
    seqdata <- seqdata[inotna,]
    seq1 <- list()
    seq2 <- list()
    if (is.null(subgroup)){
      #cat("\n No subgroup\n")
      seq1[[1]] <- seqdata[gvar==lev.g[1],]
      seq2[[1]] <- seqdata[gvar==lev.g[2],]
    }
    else {
      #cat("\n Subgroup\n")
      for (i in 1:length(lev.sg)){
        seq1[[i]] <- seqdata[gvar==lev.g[1] & sgvar==lev.sg[i],]
        seq2[[i]] <- seqdata[gvar==lev.g[2] & sgvar==lev.sg[i],]
      }
    }
  }

  # prepare samples
  G = length(seq1)
  n = matrix(NA,nrow=G,ncol=2)
  seq.a = seq.b <- list(rep(NA,G))
  for (i in 1:G) {
    if (nrow(seq1[[i]])>=nrow(seq2[[i]])) {
      n[i,1] <- nrow(seq1[[i]])
      n[i,2] <- nrow(seq2[[i]])
      seq.a[[i]] <- seq1[[i]]
      seq.b[[i]] <- seq2[[i]]
    }
    else {
      n[i,1] <- nrow(seq2[[i]])
      n[i,2] <- nrow(seq1[[i]])
      seq.a[[i]] <- seq2[[i]]
      seq.b[[i]] <- seq1[[i]]
    }
  }

  m.n = apply(n,1,max)
  n.n = apply(n,1,min)
  f.n1 <- floor(s/m.n)
  ff.n1 <- sapply(f.n1, g<-function(x){max(1,x)})
  #r.n1 = s-m.n%%s
  r.n1 = ifelse(s<m.n, s - m.n%%s, s - f.n1*m.n)
  #k.n = floor((m.n+r.n1)/n.n)
  #r.n2 = (m.n+r.n1)-k.n*n.n
  k.n = floor((ff.n1*m.n+r.n1)/n.n)
  r.n2 = (ff.n1*m.n+r.n1)-k.n*n.n

  #tmp <- data.frame(m.n,n.n,f.n1,r.n1,k.n,r.n2)
  #print(tmp)

  ## we have an error when s > min(m.n), because then we get some r.n1 > m.n
  if(any(m.n<r.n1)) {
    ii <- which(m.n<r.n1)
    print(cbind(n,m.n,r.n1))
    stop("rest r.n1 values greater than max m.n for i= ", ii, " s= ", s)
  }
  if(any(n.n<r.n2)) {
    ii <- which(n.n<r.n2)
    print(cbind(n,n.n,r.n2))
    stop("rest r.n2 values greater than min n.n for i= ", ii, " s= ", s)
  }

  r.s1=r.s2 = list(rep(NA,G))
  for (i in 1:G) {
    set.seed(seed)
    #if (s<m.n[i]) {
    #  r.s1[[i]] <- c(permute(1:m.n[i]),sample(1:m.n[i],r.n1[i],F))
    #}
    #else {
      r.s1[[i]] <- c(permute(rep(1:m.n[i],ff.n1[i])),sample(1:m.n[i],r.n1[i],F))
    #}
    #cat("length r.s1 ",length(r.s1[[i]]))
    r.s2[[i]] <- c(permute(rep(1:n.n[i],k.n[i])),sample(1:n.n[i],r.n2[i],F))
    #cat("length r.s2 ",length(r.s2[[i]]))
    r.s1[[i]] = matrix(r.s1[[i]],ncol=s)
    r.s2[[i]] = matrix(r.s2[[i]],ncol=s)
  }
  #k = rep(NA,G)
  #for (i in 1:G) k[i]<-nrow(r.s1[[i]])
  #print(k)

  nc <- ifelse(is.LRT & is.BIC, 4, 2)
  Results=matrix(NA,G,nc)

  ### new complete samples without replacement of length s over G comparisons
  for (i in 1:G) {
    t<-matrix(NA,nrow=nrow(r.s1[[i]]),ncol=nc)
    for (j in 1:nrow(r.s1[[i]])) {
      seqA<-seq.a[[i]][r.s1[[i]][j,],]
      seqB<-seq.b[[i]][r.s2[[i]][j,],]
      ##suppressMessages(t[j,]<-unlist(seq.comp2(seqA, seqB, BIC=BIC, method=method, ...)))
      suppressMessages(t[j,] <-
        seq.comp2(seqA, seqB, is.LRT=is.LRT, is.BIC=is.BIC,
        method=method, squared=squared, LRTpow=LRTpow, ...))
    }
    Results[i,]<-apply(t,2,mean)
    colnames <- NULL
    if (is.LRT) colnames <- c("LRT", "p-value")
    if (is.BIC) colnames <- c(colnames, "BIC diff.", "Bayes Factor")
    colnames(Results) <- colnames
  }
  if(!is.null(subgroup)) rownames(Results) <- lev.sg
  return(Results)
}

seq.comp2 <- function(S1,S2,is.LRT,is.BIC,method=method, squared, LRTpow,...)
{

  # compute some basic statistics
  n1 = nrow(S1)
  n2 = nrow(S2)
  S = rbind(S1,S2)
  attr(S, "weights") <- c(attr(S1, "weights"), attr(S2, "weights"))
  n0 = nrow(S)
  n.0=log((n0*n0-n0)/2)
  dist.S=dist.S1=dist.S2<-vector()


  # compute dissimilarity matrices & distances to centers
  distS <- seqdist(S, method=method, ...) #  distance matrix for overall sample
  dist.S<-disscenter(distS, squared=squared) # calculate S distance to center
  distS1 <- seqdist(S1, method=method, ...) #  distance matrix for sample 1
  distS2 <- seqdist(S2, method=method, ...) #  distance matrix for sample 2
  dist.S1<-disscenter(distS1, squared=squared) # calculate S1 distance to center
  dist.S2<-disscenter(distS2, squared=squared) # calculate S2 distance to center

  res <- NULL


  # compute LRTs and alpha probabilities
  LRT<-n0*log(sum(dist.S^LRTpow)/n0)-n0*log(sum(c(dist.S1,dist.S2)^LRTpow)/n0)

  if (is.LRT) {
    p.LRT <- pchisq(LRT,1,lower.tail=FALSE)
    res <- cbind(LRT, p.LRT)
  }

  if (is.BIC) {

    # compute BICs and adjusted BICs
    BIC <- LRT - 1*log(n0)
    # compute Bayes factors
    BF <- exp(BIC/2)

    res <- cbind(res, BIC, BF)

  # output computation results
  #results <- list("LRT comparing Groups 1 & 2"=LRT,
  #        "LRT significance level"=p.LRT,
  #        "BIC difference between Groups 1 & 2"=BIC,
  #        "Bayes factor comparing the two groups"=BF)
  }

  return(res)
}
