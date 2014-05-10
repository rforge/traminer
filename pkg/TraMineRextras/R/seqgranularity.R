seqgranularity <- function(seqdata, tspan=3, method="last"){

    metlist <- c("first","last","mostfreq")
    if (!method %in% metlist) {
        stop(" [!] method must be one of: ", paste(metlist, collapse = " "),
            call. = FALSE)
    }
    
    n <- nrow(seqdata)
    lgth <- max(seqlength(seqdata))

    new.lgth <- ceiling(lgth/tspan)
    new.lgth.f <- floor(lgth/tspan)

    newseq <- seqdata[,1:new.lgth]
    cnames <- names(seqdata)
    newcnames <- cnames[seq(from=1, to=lgth, by=tspan)]

    prev <- 0
    if (method=="first") prev <- tspan - 1
    for (i in 1:new.lgth.f){
        newseq[,i] <- seqdata[,tspan*i - prev]
    }
    if (method=="mostfreq") {
        for(i in 1:new.lgth.f) {
            print(i)
            st.freq <- seqistatd(seqdata[,(tspan*(i-1) + 1):(tspan*i)])
            alph <- colnames(st.freq)
            idmax <- apply(st.freq,1,which.max)
            st <- vector(mode="character",length(idmax))
            for (j in 1:length(idmax)){
               st[j] <- alph[idmax[j]]
            }
            newseq[,i] <- st
        }
}
    if (new.lgth > new.lgth.f){
        if (method=="first"){
            newseq[,new.lgth] <- seqdata[,tspan*new.lgth.f + 1]
        } else if (method=="mostfreq") {
            st.freq <- seqistatd(seqdata[,(tspan*new.lgth.f+1):lgth])
            alph <- colnames(st.freq)
            idmax <- apply(st.freq,1,which.max)
            st <- vector(mode="character",length(idmax))
            for (j in 1:length(idmax)){
              st[j] <- alph[idmax[j]]
            }
            newseq[,new.lgth] <- st
        } else { # method=="last"  
            newseq[,new.lgth] <- seqdata[,lgth]
        }
        newcnames[new.lgth] <- cnames[lgth]
    }
    colnames(newseq) <- newcnames
    attr(newseq,"xtstep") <- ceiling(attr(seqdata,"xtstep")/tspan)
	return(newseq)
}
