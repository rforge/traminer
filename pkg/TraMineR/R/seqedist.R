seqedist <- function(seqe, idcost, vparam, interval=TRUE, norm=TRUE){
    norm <- as.integer(norm)
    interval <- as.integer(interval)
    return(.Call(C_tmrseqedist, seqe, idcost, vparam, norm, interval));
}

seqeage <- function(seqe, eventList){
	return(.Call(C_tmreventinseq, seqe, as.integer(eventList)))
}
