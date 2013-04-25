wcKMedRange <- function(diss, kvals, weights=NULL, ...){
	if (inherits(diss, "dist")) {
		n <- attr(diss, "Size")
	}else if(is.matrix(diss)){
		n <- nrow(diss)
	} else {
		stop("[!] diss should be a squared matrix or a dist object.")
	}
	ret <- list()
	ret$kvals <- kvals
	ret$clustering <- matrix(-1, nrow=n, ncol=length(kvals))
	ret$stats <-  matrix(-1, nrow=length(kvals), ncol=10)
	i <- 1
	kendall <- .Call(wc_RClusterQualKendallFactory)
	for(k in kvals){
		cl <- wcKMedoids(diss=diss, k=k, cluster.only=TRUE, weights=weights, ...)
		ret$clustering[,i] <- cl
		stat <- wcClusterQualityInternal(diss=diss, clustering=cl, weights=weights, kendall=kendall)
		ret$stats[i,] <- stat$stats
		i <- i+1
	}
	ret$clustering <- as.data.frame(ret$clustering)
	ret$stats <- as.data.frame(ret$stats)
	colnames(ret$stats) <- names(stat$stats)
	colnames(ret$clustering) <- paste("cluster", ret$kvals, sep="")
	rownames(ret$stats) <- paste("cluster", ret$kvals, sep="")
	class(ret) <- c("clustrange", class(ret))
	return (ret)
}

as.clustrange <- function(object, diss, weights=NULL, ...){
	UseMethod("as.clustrange")
}

as.clustrange.hclust <- function(object, diss, weights=NULL, ncluster, ...){
	
	if(ncluster<3){
		stop(" [!] ncluster should be greater than 2.")
	}
	if (is.null(n1 <- nrow(object$merge)) || n1 < 1) {
		stop("invalid 'object' (merge component)")
	}
    n <- n1 + 1
	if(ncluster > n){
		stop(" [!] ncluster should be less than ", n)
	}
	
	pred <- data.frame(Split2=factor(cutree(object, 2)))
	for(p in 3:ncluster){
		pred[, paste("Split", p, sep="")] <- factor(cutree(object, p))
	}
	object <- pred
	as.clustrange(object, diss=diss, weights=weights, ...)
}

as.clustrange.twins <- function(object, diss, weights=NULL, ncluster, ...) {
	return(as.clustrange.hclust(object, diss=diss, weights=weights, ncluster=ncluster, ...))
}

as.clustrange.default <- function(object, diss, weights=NULL, ...){
	ret <- list()
	ret$clustering <- as.data.frame(object)
	numclust <- ncol(ret$clustering)
	ret$kvals <- numeric(numclust)
	ret$stats <-  matrix(-1, nrow=numclust, ncol=10)
	## print("BuildingKendall")
	kendall <- .Call(wc_RClusterQualKendallFactory)
	## print(class(kendall))
	## print(kendall)
	## print("Kendall")
	for(i in 1:numclust){
		ret$kvals[i] <- length(unique(ret$clustering[,i]))
		## print("Starting")
		cl <- wcClusterQualityInternal(diss, ret$clustering[,i], weights=weights, kendall=kendall)
		ret$stats[i,] <- cl$stats
	}
	ret$stats <- as.data.frame(ret$stats)
	colnames(ret$stats) <- names(cl$stats)
	colnames(ret$clustering) <- paste("cluster", ret$kvals, sep="")
	rownames(ret$stats) <- paste("cluster", ret$kvals, sep="")
	class(ret) <- c("clustrange", class(ret))
	return(ret)
}
print.clustrange <- function(x, digits=2, ...){
	x <- round(x$stats, digits)
	print(x, ...)
}

normalize.values.all <- function(stats, norm){
	for(i in 1:ncol(stats)){
		stats[,i] <- normalize.values(stats[, i], norm)
	}
	return(stats)
}
normalize.values <- function(stats, norm){
	if(norm == "range") {
		stats <- (stats-min(stats))/(max(stats)-min(stats))
	} else if (norm=="zscore") {
		stats <- (stats - mean(stats ))/(sqrt(var(stats)))
	}else if (norm=="zscoremed") {
		stats <- (stats - median(stats))/(median(abs(stats - median(stats))))
	}
	return(stats)
}


plot.clustrange <- function(x, stat="noCH", legendpos="bottomright", norm="none", withlegend=TRUE, lwd=1, col=NULL, ylab="Indicators", xlab="N clusters", ...){
	kvals <- x$kvals
	if(length(stat)==1){
		if(stat=="all"){
			stats <- x$stats
			if(norm=="none"){
				norm <- "range"
			}
		}else if(stat=="noCH"){
			stats <- x$stats[,-grep("CH", colnames(x$stats))]
		}
		else{
			stats <- x$stats[,stat]
		}
	}else{
		if("RHC" %in% stat){
			stats <- x$stats[ , stat[stat != "RHC"]]
			stats[, "RHC"] <- 1 - x$stats[ , "HC"]
		}
		else{
			stats <- x$stats[, stat]
		}
	}
	stats <- normalize.values.all(stats, norm)
	ylim <- range(unlist(stats), finite=TRUE)
	plot(kvals, xlim=range(kvals, finite=TRUE), ylim=ylim, type="n", ylab=ylab, xlab=xlab, ...)
	labels <- character(ncol(stats))
	if(is.null(col)){
		allnames <- colnames(x$stats)
		cols <- brewer.pal(length(allnames)+1, "Set3")[-2]
		names(cols) <- allnames
		cols["RHC"] <- cols["HC"]
		cols <- cols[colnames(stats)]
	} else {
		if(length(col) < ncol(stats)){
			stop(" [!] You should specify at least one color per quality measure to plot.")
		}
		cols <- col
	}
	for(i in 1:ncol(stats)){
		ss <- stats[,i]
		lines(kvals, ss, col=cols[i], lwd=lwd, ...)
		labels[i] <- paste(colnames(stats)[i], "(", round(min(stats[,i]), 2),"/", round(max(stats[,i]),2), ")")
	}
	if(withlegend) {
		legend(legendpos, fill=cols[1:ncol(stats)], legend=labels)
	}
}

summary.clustrange <- function(object, max.rank=1, ...){
	nstat <- ncol(object$stats)
	clusterrank <- matrix(NA, nrow=nstat, ncol=max.rank*2)
	rownames(clusterrank) <- colnames(object$stats)
	clindices <- ((1:max.rank) *2) -1
	valindices <- ((1:max.rank) *2)
	cnames <- character(max.rank*2)
	cnames[clindices] <- paste(1:max.rank, "N groups", sep=". ")
	cnames[valindices] <- paste(1:max.rank, " stat", sep=". ")
	colnames(clusterrank) <- cnames
	for(s in colnames(object$stats)){
		od <- order(object$stats[, s], decreasing=(s!="HC"))[1:max.rank]
		clusterrank[s, clindices] <- object$kvals[od]
		clusterrank[s, valindices] <- object$stats[od, s]
	}
	
	return(as.data.frame(clusterrank, check.names=FALSE))
}
