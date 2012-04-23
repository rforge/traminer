seqgranularity <- function(seqdata, tspan=3, method="last"){

    n <- nrow(seqdata)
    lgth <- max(seqlength(seqdata))

    new.lgth <- ceiling(lgth/tspan)
    new.lgth.f <- floor(lgth/tspan)

    newseq <- seqdata[,1:new.lgth]
    cnames <- names(seqdata)
    newcnames <- cnames[seq(from=1, to=lgth, by=tspan)]
    for (i in 1:new.lgth.f){
        newseq[,i] <- seqdata[,tspan*i]
    }
    if (new.lgth > new.lgth.f){
        newseq[,new.lgth] <- seqdata[,lgth]
        newcnames[new.lgth] <- cnames[lgth]
    }
    colnames(newseq) <- newcnames
    attr(newseq,"xtstep") <- ceiling(attr(seqdata,"xtstep")/tspan)
	return(newseq)
}

