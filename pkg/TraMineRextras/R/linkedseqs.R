### linked sequences function for linked polyadic sequence data
### Based on Version 1.0.1 (14.12.20), Tim Liao, University of Illinois
### the function uses sequence data defined by TraMinerR
### Implemented in TraMineRextras by Gilbert Ritschard

linkedseqs <- function (seqlist, a=1, method="HAM", ..., w=rep(1,ncol(combn(1:length(seqlist),2))),
              s=36963, T=1000, core=1) {
  #require(doParallel)
  #require(TraMineR)
  cl <- makeCluster(core, type="SOCK")
  registerDoParallel(cl)

  set.seed(s)
  P <- length(seqlist)
  c <- combn(1:P,2)
  d = ncol(c)
  n = nrow(seqlist[[1]])
  s = ncol(seqlist[[1]])
  p.seq=p.temp=p1.temp <- list()
  polyads.dist=polyads.dummy=test.p <- array(0,c(1,n))

  for (p in 1:P) p.seq[[p]] <- suppressMessages(seqformat(seqlist[[p]],to="STS"))

  l <- rep(NA,P)
  random.dist <- foreach (i=1:T, .combine='c', .packages='TraMineR') %dopar% {
    for (p in 1:P) {
      l[p] <- sample(1:nrow(seqlist[[p]]),1)
      if (a==1) p1.temp[[p]] <- seqlist[[p]][l[p],]
      if (a==2) {
        p.seq[[p]] <- suppressMessages(seqformat(seqlist[[p]],to="STS"))
        p1.temp[[p]] <- seqdef(sample(p.seq[[p]][l[p],],s))
      }
    }

    seq.temp <- do.call("rbind",p1.temp)
    sum(suppressMessages(seqdist(seq.temp,method=method,...))[upper.tri(matrix(NA,P,P))]*w)/d
  }
  stopCluster(cl)

  for (j in 1:n) {
    for (p in 1:P) p.temp[[p]] <- seqlist[[p]][j,]
    seq.temp <- suppressMessages(seqdef(do.call("rbind",p.temp)))
    polyads.dist[j] <- sum(suppressMessages(seqdist(seq.temp,method=method,...))[upper.tri(matrix(NA,P,P))]*w)/d
    test.p[j] <- sum(polyads.dist[j]<random.dist)/T
    if (test.p[j]>0.95) polyads.dummy[j] <- 1
    else polyads.dummy[j] <- 0
  }
  mean.U <- mean(random.dist) - polyads.dist
  U.p <- 2*pt(mean.U/{sd(random.dist)/sqrt(T)},T-1,lower.tail=F)
  list(mean.obs=mean(polyads.dist),U=mean.U,U.tp=U.p,V=test.p,
       V.95=polyads.dummy,observed.dist=polyads.dist,random.dist=random.dist)
}
