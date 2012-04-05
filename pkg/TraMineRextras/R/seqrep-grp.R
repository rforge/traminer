seqrep.grp <- function(seqdata, group=NULL, mdis=NULL, ret="stat", ...){

    if (!(ret %in% c("stat","rep","both"))){
        stop("ret should be one of 'stat', 'rep' or 'both'",
            call. = FALSE)
        }
    grp <- group
    if (is.null(grp)) grp <- rep(1, nrow(seqdata))
    if (length(grp) != nrow(seqdata)){
        stop("length(grp) not equal to number of sequences",
            call. = FALSE)
        }

    levg <- levels(grp <- factor(grp))

    if (is.null(mdis)) mdis <- seqdist(seqdata, method="LCS")
    dmax <- max(mdis)

    q.gr <- gr <- list()
    for (i in 1:length(levg))
        {
        ig <- which(grp==levg[i]) 
        gr[[i]] <- seqrep(seqdata[ig,], dist.matrix=mdis[ig,ig], dmax=dmax, ...)
        q.gr[[i]] <- attr(gr[[i]],"Statistics")
        }
    names(q.gr) <- names(gr) <- levg
    if (ret=="stat") return(q.gr)
    if (ret == "rep") return(gr)
    if (ret=="both") return(list(gr, q.gr))
}

