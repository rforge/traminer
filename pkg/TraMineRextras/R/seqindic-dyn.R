## Evolution of indicator
## Based on Pelletier et al., 2020

seqindic.dyn <- function(seqdata, indic="cplx", window.size = .2,
      with.missing=TRUE, ...) {

	if (!inherits(seqdata,"stslist"))
		TraMineR:::msg.stop("data is NOT a sequence object, see seqdef function to create one")

  if (length(indic) > 1) {
    TraMineR:::msg.warn("Vector indic, only first value is used!")
    indic <- indic[1]
  }
  group.list <- c("all","basic","diversity","complexity","binary","ranked")
  if (indic %in% group.list)
		TraMineR:::msg.stop("Bad indic value, group name not supported!")

  lgth <- seqlength(seqdata, with.missing=TRUE)
  maxl <- max(lgth)
  if (window.size <= 0)
    TraMineR:::msg.stop("bad window.size value!")
  if (window.size < 1)
    step <- ceiling(window.size * maxl)
  else {
    step <- min(window.size, maxl)
  }

  re <- step
  rs <- 1

  nr <- nrow(seqdata)
  nc <- maxl-step+1
  ind.dyn <- matrix(NA, nrow=nr, ncol=nc)
  colnames(ind.dyn) <- colnames(seqdata)[ceiling(step/2):(nc+ceiling(step/2)-1)]
  rownames(ind.dyn) <- rownames(seqdata)
  j <- 1

  while (re < maxl + 1) {
    ind.dyn[,j] <- seqindic(seqdata[,rs:re], indic=indic, with.missing=with.missing, ...)[,1]
    j  <- j+1
    rs <- rs+1
    re <- re+1
  }


  class(ind.dyn) <- c("dynin",class(ind.dyn))
  attr(ind.dyn,"xtstep") <- attr(seqdata,"xtstep")
  attr(ind.dyn,"tick.last") <- attr(seqdata,"tick.last")
  attr(ind.dyn, "window.size") <- step
  attr(ind.dyn, "indic") <- indic

  return(ind.dyn)
}

#################################

print.dynin <- function(x, ...){
  dimx <- dim(x)
  dimnamesx <- dimnames(x)
  attributes(x) <- NULL
  dim(x) <- dimx
  dimnames(x) <- dimnamesx
  print(x, ...)
}

#################################

plot.dynin <- function(x, fstat=mean, group=NULL,
     main=NULL, col=NULL, lty=NULL, lwd=3.5, ylim=NULL,
     ylab=NULL, xlab=NULL, xtlab=NULL, xtstep=NULL, tick.last=NULL,
     with.legend=TRUE, glabels=NULL, legend.pos="topright",
     horiz=FALSE, cex.legend=1, ret.table=FALSE, ...){

  if (class(fstat)!="function") TraMineR:::msg.stop("stat must be a function!")

  nc <- ncol(x)
  nr <- nrow(x)
  if (is.null(group)){
    group <- rep(1,nr)
  }
  if (!is.factor(group)) group <- factor(group)
  ngrp <- length(levels(group))

  tab.grp <- matrix(NA, nrow=ngrp, ncol=nc)
  colnames(tab.grp) <- colnames(x)
  rownames(tab.grp) <- levels(group)

  #for(i in 1:nc) {
  #  tab.grp[,i] <- tapply(x[,i],group,fstat)
  #}

  tab.grp <- apply(x,2,tapply,group,fstat)

  #default.col <- brewer.pal(9,"Set1")
  default.col <- qualitative_hcl(ngrp, palette = "Dark 3")
  ##default.col <- c("red","blue","black","magenta","green")

  if(is.null(col)) {
     #col <- colors.list[1:k]
     col <- default.col
  }
  kk <- ceiling(ngrp/length(col))
  col <- rep(col,kk)
  col <- col[1:ngrp]

  default.lty <- c("solid","dashed","dotted")
  if(is.null(lty)) {
     lty <- default.lty
  }
  kk <- ceiling(ngrp/length(lty))
  lty <- rep(lty,kk)
  lty <- lty[1:ngrp]

  kk <- ceiling(ngrp/length(lwd))
  lwd <- rep(lwd,kk)
  lwd <- lwd[1:ngrp]

  if(is.null(glabels)) glabels <- levels(group)
  if(is.null(xtstep)) xtstep <- attr(x,"xtstep")
  if(is.null(tick.last)) tick.last <- attr(x,"tick.last")
  if(is.null(xtlab)) xtlab <- colnames(x)
  if(is.null(xlab)) xlab<-"Window center"

  if(is.null(ylab)) ylab<-attr(x,"indic")
  main=paste("Dynamic index, window =",attr(x,"window.size"))

  if(is.null(ylim)){
        maxe <- max(tab.grp, na.rm=TRUE)
        mine <- min(tab.grp, na.rm=TRUE)
        if (mine==maxe) maxe <- mine + .1
        ylim <- c(floor(10*mine),ceiling(10*maxe))/10
  }


  plot(0, type= "n", axes=FALSE, xlab=xlab, ylab=ylab, main=main, ylim=ylim, xlim=c(1,nc), ...)
  for (i in 1:ngrp) {
     lines(tab.grp[i,], col=col[i],  type="l", lty=lty[i], lwd=lwd[i], ...)
  }
	tpos <- seq(from=1, to=nc, by=xtstep)
  if (tick.last & tpos[length(tpos)] < nc) tpos <- c(tpos,nc)
  axis(1,labels=xtlab[tpos],at=tpos)
  #axis(1)
  axis(2)
  non.na.rows <- apply(tab.grp,1,function(x)all(!is.na(x)))
  if(with.legend & ngrp>1){
    legend(legend.pos, legend=glabels[non.na.rows],
      lwd=lwd[non.na.rows], lty=lty[non.na.rows], col=col[non.na.rows],
      horiz=horiz, cex=cex.legend)
  }

  if (ret.table) return(tab.grp)
}
