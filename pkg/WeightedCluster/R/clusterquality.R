wcClusterQuality <- function(diss, clustering, weights=NULL){
	if (inherits(diss, "dist")) {
		isdist <- TRUE
		nelements <- attr(diss, "Size")
	}else if(is.matrix(diss)){
		isdist <- FALSE
		nelements <- nrow(diss)
		if(ncol(diss)!=nelements){
			stop("[!] diss should be a squared matrix or a dist object.")
		}
	} else {
		stop("[!] diss should be a squared matrix or a dist object.")
	}
	if(is.null(weights)){
		weights <- as.double(rep(1.0, nelements))
	}
	clusterF <- factor(clustering)
	if(nlevels(clusterF) < 2){
		stop("[!] The clustering should have at least two different values.")
	}
	clustering <- as.integer(as.integer(clusterF)-1)
	if(length(clustering)!=nelements|| length(weights)!=nelements){
		stop("[!] different number of elements in diss, clustering and/or weights arguments.")
	}
	
	ncluster <- max(clustering)+1
	cq <- .Call("RClusterQual", diss, clustering, as.double(weights), as.integer(ncluster), as.integer(isdist), as.integer(0))
	names(cq[[1]]) <-c("PBC", "HG", "HGSD", "ASW", "CH", "R2", "CHsq", "R2sq", "HC")
	names(cq[[2]]) <-levels(clusterF)
	names(cq) <- c("stats", "ASW")
	if(any(xtabs(weights~clustering)<1)){
		cq$stats
		warning(" [!] silhouette can not be computed for clusters with less than one observation. Use silhouette.weight='weight'.\n")
	}
	return(cq)
	#define ClusterQualHPG 0
#define ClusterQualHG 1
#define ClusterQualHGSD 2
#define ClusterQualASW 3
#define ClusterQualF 4
#define ClusterQualR 5
#define ClusterQualF2 6
#define ClusterQualR2 7
}

wcSilhouetteObs <- function(diss, clustering, weights=NULL, silhouette.weight="replicate"){
	if (inherits(diss, "dist")) {
		isdist <- TRUE
		nelements <- attr(diss, "Size")
	}else if(is.matrix(diss)){
		isdist <- FALSE
		nelements <- nrow(diss)
		if(ncol(diss)!=nelements){
			stop("[!] diss should be a squared matrix or a dist object.")
		}
	} else {
		stop("[!] diss should be a squared matrix or a dist object.")
	}
	if(is.null(weights)){
		weights <- as.double(rep(1.0, nelements))
	}
	clusterF <- factor(clustering)
	if(nlevels(clusterF) < 2){
		stop("[!] The clustering should have at least two different values.")
	}
	clustering <- as.integer(as.integer(clusterF)-1)
	if(length(clustering)!=nelements|| length(weights)!=nelements){
		stop("[!] different number of elements in diss, clustering and/or weights arguments.")
	}
	if(silhouette.weight != "replicate"){
		if(any(xtabs(weights~clustering)<1)){
			stop(" [!] silhouette can not be computed for clusters with less than one observation. Use silhouette.weight='weight'.\n")
		}
	}
	silhouette.weight <- as.integer(silhouette.weight !="replicate")
	
	ncluster <- max(clustering)+1
	asw <- .Call("RClusterComputeIndivASW", diss, clustering, as.double(weights), as.integer(ncluster), as.integer(isdist), as.integer(silhouette.weight))
	return(asw)
	#define ClusterQualHPG 0
#define ClusterQualHG 1
#define ClusterQualHGSD 2
#define ClusterQualASW 3
#define ClusterQualF 4
#define ClusterQualR 5
#define ClusterQualF2 6
#define ClusterQualR2 7
}


clusterqualitySimpleBoot <- function(diss, clustering, weights=NULL, R=999, replicates=TRUE, ...){
	if (inherits(diss, "dist")) {
		isdist <- TRUE
		nelements <- attr(diss, "Size")
	}else if(is.matrix(diss)){
		isdist <- FALSE
		nelements <- nrow(diss)
		if(ncol(diss)!=nelements){
			stop("[!] diss should be a squared matrix or a dist object.")
		}
	} else {
		stop("[!] diss should be a squared matrix or a dist object.")
	}
	if(is.null(weights)){
		weights <- as.double(rep(1.0, nelements))
	}
	clusterF <- factor(clustering)
	if(nlevels(clusterF) < 2){
		stop("[!] The clustering should have at least two different values.")
	}
	clustering <- as.integer(as.integer(clusterF)-1)
	if(length(clustering)!=nelements|| length(weights)!=nelements){
		stop("[!] different number of elements in diss, clustering and/or weights arguments.")
	}
	ncluster <- max(clustering)+1
	totweights <- sum(weights)
	fenv <- new.env()
	fenv$R <- 0
	if(replicates){
		prob <- weights/totweights
		sampleint <- as.integer(floor(totweights))
		internalsample <- function(){
			return(as.integer(sample.int(nelements, size=sampleint, replace=TRUE, prob=prob)-1L))
		}
		t <- .Call("RClusterQualSimpleBoot", diss, clustering, as.double(weights), as.integer(ncluster), as.integer(R),  quote(internalsample()), environment(), sampleint)
		colnames(t) <-c("PBC", "CH", "R2", "CHsq", "R2sq")
		t0 <- t[1,]
		bts <- list(t0=t0, t=t)
	} else {
		SimpleBootStatistic <- function(d, w) {
			if(fenv$R==0){
				w <- as.double(weights)
				fenv$R <- 1
			}else{
				w <- as.double(w*totweights)
			}
			cq <- .Call("RClusterQualSimple", diss, clustering, w, as.integer(ncluster))
			names(cq) <-c("PBC", "HG", "HGSD", "ASW", "CH", "R2", "CHsq", "R2sq", "HC")
			return(cq[c(1, 5:8)])
		}
		bts <- boot(clustering, stype="w", statistic=SimpleBootStatistic, weights=weights, R=R, ...)
	}
	return(bts)
}