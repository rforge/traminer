## ====================================================
## Computing distances between sequences
## ====================================================

seqdistOO <- function(seqdata, method, refseq=NULL, norm=FALSE,
        indel=1, sm=NA, with.missing=FALSE, full.matrix=TRUE,
        kweights=rep(1.0, ncol(seqdata)), tpow=1,
        expcost=0.5, context=1-2*expcost,
        link="mean", h=0.5, nu=0,
        transindel="constant", otto, previous=FALSE,
        addcolumn=TRUE, numcpu=1,
        breaks=NULL, step=1, overlap=FALSE, weighted=TRUE) {
	gc(FALSE)
	warnings(" [>] seqdistOO is still under development. Some configuration of parameters may produce errors.")
    ## Checking method name
	metlist <- c("OM", "OMopt", "LCP", "LCS", "LCSopt", "RLCP", "DHD", "HAM",
                 "NMS", "NMSMST", "SVRspell",
                 "OMloc", "OMspell", "OMslen", "OMstran", "TWED")
	OMvariant <- c("OM", "OMloc", "OMspell", "OMslen")
	if (missing(method)) {
		stop(" [!] You should specify a method to compute the distances. It must be one of: ", paste(metlist,collapse=" "))
	}
	if (!method %in% metlist) {
		stop(" [!] Method must be one of: ", paste(metlist, collapse=" "), call.=FALSE)
	}
	debut <- proc.time()
	if(method=="OM" && length(indel)>1){
		method <- "OMloc"
	}
	## Checking correct arguments
	if (!inherits(seqdata,"stslist")) {
		stop(" [!] data is not a state sequence object, use 'seqdef' function to create one", call.=FALSE)
	}
	if (method=="SVRspell") {
		method <- "NMSMST"
		SVRspell <- TRUE
	} else{
		SVRspell <- FALSE
	}
	if (method=="OMopt") {
		method <- "OM"
		optimized <- TRUE
	}
	else if (method=="LCSopt") {
		method <- "LCS"
		optimized <- TRUE
	} else if(method %in% c("CHI2", "EUCLID", "OMstran")){
		if(method=="CHI2"){
			distances <- chisqdist(seqdata, breaks=breaks, step=step, with.missing=with.missing, norm=norm, weighted=weighted, overlap=overlap)
		} else if(method=="EUCLID"){
			distances <- chisqdist(seqdata, breaks=breaks, step=step, with.missing=with.missing, norm=norm, weighted=FALSE, overlap=overlap, euclid=TRUE)
		} else if(method=="OMstran"){
			distances <- seqdistOMTrans(seqObj=seqdata, indel=indel, sm=sm, transmetric=transindel, stateweight=otto, with.missing=with.missing, norm=norm, previous=previous, addcolumn=addcolumn, full.matrix=full.matrix, weighted=weighted, numcpu=numcpu)
		}
		if (full.matrix && inherits(distances, "dist")) {
			return(TraMineR:::dist2matrix(distances))
		}
		else {
			return(distances)
		}
	}
	else {
		optimized <- FALSE
	}
	## Taking care of correct normalization settings
	if (is.logical(norm)) {
		if (norm) {
		## Normalize using Elzinga for LCP, LCS, RLCP
			if (method %in% c("LCP", "LCS", "RLCP")) {
				norm <- 2
			} else {## Normalize using Abbott for OM, HAM, and DHD
				norm <- 1
			}
      	} else {
      		norm <- 0
		}
	} else if (is.character(norm)) {
		## Using normalization name
		## Cast to integer for c code and normdist function
		## Match return the position, removing 1 to start at zero
		normIndice <- match(norm, c("none", "maxlength", "gmean", "maxdist", "YujianBo")) -1

		if (is.na(normIndice)) {  ##Not found
			stop(" [!] Unknown distance normalization method ", norm)
		}
		norm <- normIndice
	} else {
		stop(" [!] Unknown distance normalization method ", norm)
	}
	## Checking missing values
	if (!with.missing && any(seqdata==attr(seqdata,"nr"))) {
		stop("found missing values in sequences, please set 'with.missing=TRUE' to nevertheless compute distances")
	}
	if (method == "OM" && is.character(sm)) {
		if (sm == "TRATE" || sm =="CONSTANT") {
			sm <- seqsubm(seqdata, method=sm, cval=2, with.missing=with.missing, miss.cost=2)
		} else {
			stop(" [!] Unknown method ", sm, " to compute substitution costs")
		}
	}
	
	## =====================
	## Base information
	## =====================
	
	n <- nrow(seqdata)
	alphabet <- attr(seqdata,"alphabet")
	alphsize <- length(alphabet)
	message(" [>] ",n," sequences with ", alphsize, " distinct events/states")
	## Gaps in sequences
	if (with.missing) {
		alphabet <- c(alphabet,attr(seqdata,"nr"))
		alphsize <- length(alphabet)
		message(" [>] including missing value as additional state" )
	}
	
	## Checking methods that are treated the same
	methodname <- method
	if (method == "LCS") {
		method <- "OM"
		sm <- suppressMessages(seqsubm(seqdata, method="CONSTANT", cval=2, with.missing=with.missing, miss.cost=2))
		indel <- 1
	} else if (method == "HAM") {
		method <- "DHD"
		
		if (!is.null(dim(sm))) {
			TraMineR.checkcost(sma=sm, seqdata=seqdata, with.missing=with.missing)
			if(is.matrix(sm)){
				costs <- array(0, dim=c(alphsize, alphsize, ncol(seqdata)))
				for(i in 1:ncol(seqdata)){
					costs[,,i] <- sm
				}
				sm <- costs
			}
		} else {
			sm <- suppressMessages(seqsubm(seqdata, "CONSTANT", cval=1, with.missing=with.missing,
				miss.cost=1, time.varying=TRUE))
		}
	}


	## ===========================
	## Checking correct size of sm
	## ===========================
	
	## Checking if substitution cost matrix contains values for each state
	## and if the triangle inequality is respected
	if (method %in% OMvariant) {
		TraMineR.checkcost(sma=sm, seqdata=seqdata, with.missing=with.missing, indel=indel)
	}
	## Checking if substitution cost matrix contains values for each state
	## and if the triangle inequality is respected
	if(methodname == "DHD"){
		## User entered substitution cost
		if(length(sm)>1){
			TraMineR.checkcost(sma=sm, seqdata=seqdata, with.missing=with.missing)
		}
		else {
			sm <- seqsubm(seqdata, "TRATE", cval=4, with.missing=with.missing,
					miss.cost=4, time.varying=TRUE)
		}
	}


	## ==============
	## Preparing data
	## ==============
	if(method %in% c("NMSMST", "OMspell")) {
		seqduree <- seqdur(seqdata, with.missing=with.missing)^tpow
		
		seqdatadss <- seqdss(seqdata, with.missing=with.missing)
	}else{
		seqduree <- NULL
	}
	seqdata <- seqnum(seqdata, with.missing=with.missing)

	## Selecting distinct sequences only and saving the indexes
	dseq <- unique(seqdata)
	mcorr <- match(seqconc(seqdata), seqconc(dseq))
	
	if(method %in% c("NMSMST", "OMspell")) {
		mcorrReverse <- match(seqconc(dseq), seqconc(seqdata))
		seqduree <- seqduree[mcorrReverse,]
		dseq <- seqnum(seqdatadss[mcorrReverse,], with.missing=with.missing)
	}
	nd <- nrow(dseq)
	message(" [>] ", nd," distinct sequences")
	
	slength <- seqlength(dseq)
	if (method=="OMslen") {
		seqduree <- matrix(0, nrow=nrow(dseq), ncol=ncol(dseq))
		ddur <- seqdur(dseq)
		for(i in 1:nrow(dseq)){
			y <- ddur[i, !is.na(ddur[i,])]
			seqduree[i, 1:sum(y)] <- rep(y, times=y)
		}
		seqduree <- seqduree^(-1*h)
	}

	dseq <- TraMineR:::seqasnum(dseq, with.missing=with.missing)

	message(" [>] min/max sequence length: ",min(slength),"/",max(slength))

	## ===================================
	## Preparing data for Hamming distance
	## ===================================
	if (method=="DHD") {
		if(length(unique(slength))>1) {
			## Hamming is not defined for sequence of different length
			stop(methodname, " distance can only be computed between sequences of equal length")
		}
		## Here we use the indel parameter for the C function
		## But it contain the maximum possible cost of the hamming distance
		indel <- 0
		for (i in 1:max(slength)) {
			indel <- indel + max(sm[,,i])
		}
	}
	message(" [>] computing distances using ", methodname,
		ifelse(norm!=0," normalized", ""), " metric")
	
	
	params <- list()
	if (method=="OM") {
		## One for OM, 2 for LCP
   		disttype <- as.integer(1)
		params[["scost"]] <- sm
		params[["alphasize"]] <- alphsize
		params[["indel"]] <- indel
	}
	else if (method=="LCP") {
		disttype <- as.integer(2) ## One for OM, 2 for LCP
		params[["sign"]] <- as.integer(1)
	}
	else if (method=="RLCP") {
		disttype <- as.integer(3) ## One for OM, 2 for LCP
		params[["sign"]] <- as.integer(-1)
	}
	else if (method=="DHD") {
		disttype <- as.integer(4) ## 4 for DHD or HAM
		params[["scost"]] <- sm
		params[["alphasize"]] <- alphsize
		params[["maxdist"]] <- indel
	}
	else if (method=="NMS") {
		## One for OM, 2 for LCP
   		disttype <- as.integer(5)
		kweights2 <- double(ncol(dseq))
		kweights2[1:min(ncol(dseq),length(kweights))] <- kweights[1:min(ncol(dseq),length(kweights))]
		params[["kweights"]] <- as.double(kweights2)
		params[["distMethod"]] <- as.integer(2)
		params[["distTransform"]] <- as.integer(link=="log")
		if(!is.na(sm)){
			disttype <- as.integer(12)
			params[["softmatch"]] <- sm
			params[["alphasize"]] <- alphsize
			if (nrow(sm)!=alphsize | ncol(sm)!=alphsize) {
				stop(" [!] size of soft matching matrix must be ", alphsize,"x", alphsize)
			}
			eg <- eigen(sm)
			if(any(eg$values<0)){
				warning(" [!] matching between states is not semipositive definite. Eigen values: ", paste(eg$values, collapse=" "))
			}
		}
	}
	else if (method=="NMSMST") {
		## One for OM, 2 for LCP
		kweights2 <- double(ncol(dseq))
		kweights2[1:min(ncol(dseq),length(kweights))] <- kweights[1:min(ncol(dseq),length(kweights))]
   		disttype <- as.integer(6)
		params[["seqdur"]] <- as.double(seqduree)
		params[["kweights"]] <- as.double(kweights2)
		params[["distMethod"]] <- as.integer(2)
		params[["distTransform"]] <- as.integer(link=="log")
		if(length(sm)>1||SVRspell){
			if(SVRspell){ ## SVRspell
				disttype <- as.integer(13)
				if(length(sm)==1){
					sm <- matrix(0, ncol=alphsize, nrow=alphsize)
					diag(sm) <- 1
				}
			} else {
				disttype <- as.integer(11)
			}
			params[["softmatch"]] <- sm
			params[["alphasize"]] <- alphsize
			if (nrow(sm)!=alphsize | ncol(sm)!=alphsize) {
				stop(" [!] size of soft matching matrix must be ", alphsize,"x", alphsize)
			}
			eg <- eigen(sm)
			if(any(eg$values < -1e-07)){
				warning(" [!] matching between states is not semipositive definite. Eigen values: ", paste(round(eg$values, 3), collapse=" "))
			}
		}
	}else if (method=="OMloc") {
		## One for OM, 2 for LCP
   		disttype <- as.integer(7)
		params[["scost"]] <- sm
		params[["alphasize"]] <- alphsize
		params[["indel"]] <- indel
		if(length(indel)>1){
			params[["indels"]] <- indel
			params[["indelmethod"]] <- as.integer(0)
		}else  {
			params[["indels"]] <- rep(indel, alphsize)
			if(link=="min"){
				params[["indelmethod"]] <- as.integer(2)
			} else if(link=="mean"){
				## Maximum indel cost
				params[["indel"]] <- max(sm)*(expcost+context)
				params[["indelmethod"]] <- as.integer(1)
			} else if(link=="previous"){
				disttype <- as.integer(9)
				params[["firststatemethod"]] <- as.integer(1)
			}
			else if(link=="previous_nofirst"){
				disttype <- as.integer(9)
				params[["firststatemethod"]] <- as.integer(0)
			}
			
			message(" [>] 2*expcost+context=", 2*expcost+context)
		}
		params[["localcost"]] <- context
		params[["timecost"]] <- expcost
		
	}
	else if (method=="OMspell") {
		## One for OM, 2 for LCP
		params[["seqdur"]] <- as.double(seqduree-1)
		disttype <- as.integer(8)
		params[["scost"]] <- sm
		params[["alphasize"]] <- alphsize
		params[["indel"]] <- indel
		params[["timecost"]] <- expcost
	}
	else if (method=="OMslen") {
		## One for OM, 2 for LCP
		params[["seqdur"]] <- as.double(seqduree)
		disttype <- as.integer(10)
		if(link=="mean"){
			params[["sublink"]] <- as.integer(1)
			params[["scost"]] <- sm/2
		}
		else if (link=="gmean"){
			params[["sublink"]] <- as.integer(0)
			params[["scost"]] <- sm
		}
		else if (link=="max"){
			params[["sublink"]] <- as.integer(2)
			params[["scost"]] <- sm
		}
		params[["alphasize"]] <- alphsize
		params[["indel"]] <- indel
	} else if (method=="TWED") {
		## One for OM, 2 for LCP
   		disttype <- as.integer(14)
		params[["scost"]] <- sm
		params[["alphasize"]] <- alphsize
		params[["indel"]] <- indel
		params[["lambda"]] <- h
		params[["nu"]] <- nu
	}else {
		stop(" [!] unsupported distance method ", method)
	}
	
	## Function and arguments
	if (!missing(refseq) && !is.null(refseq)) {
		if (is.numeric(refseq) && refseq>0 && refseq <= nrow(seqdata)) {
			message(" [>] using sequence ", refseq,": ",
				suppressMessages(seqformat(seqdata[refseq,], from="STS", to="SPS", compressed=TRUE)), " as reference")
			refseqid <- mcorr[refseq]
		} else {
			stop("[!] invalid reference sequence", call.=FALSE)
		}
		## Getting refseq
		##User specified
		# if (inherits(refseq,"stslist") && nrow(refseq)==1) {
			# compseq <- refseq
			# message(" [>] using (external) sequence ",
				# suppressMessages(seqformat(compseq, from="STS", to="SPS", compressed=TRUE)), " as reference")
		# }
		# ## Most frequent sequence as reference
		# else if (refseq==0) {
			# mfseq <- seqtab(seqdata, tlim=1)
			# message(" [>] using most frequent sequence as reference: ",
				# suppressMessages(seqformat(mfseq, from="STS", to="SPS", compressed=TRUE)))
			# idxmfseq <- suppressMessages(seqfind(mfseq, seqdata))
			# message(" [>] most frequent sequence appears ", length(idxmfseq), " times")
			# compseq <- seqdata[idxmfseq[1],]
		# }
		# ## Indice of sequence given as reference
		# else
		# ## Length of compseq
		## lcompseq <- seqlength(compseq)
		## Vector of distance
		## m <- vector(mode="numeric", length=nd)
		## compseq <- TraMineR:::seqasnum(seqnum(compseq), with.missing=with.missing)
		## ddseq <- rbind(dseq[1,], dseq)
		## ddseq[1,1:lcompseq] <- compseq[1:lcompseq]
		## slength <- c(lcompseq, slength)
		## SEXP cstringrefseqdistanceOO(SEXP Ssequences, SEXP seqdim, SEXP lenS, SEXP paramS, SEXP normS, SEXP disttypeS, SEXP refseqS) {
		distances <- .Call(TMR_cstringrefseqdistanceOO,
			as.integer(dseq),
			as.integer(dim(dseq)),
			as.integer(slength),
			params,
			as.integer(norm),
			disttype,
			as.integer(refseqid))
		## distances <- distances[2:length(distances)]
		## Constructing the final distance vector
		#mcorr <- match(seqconc(seqdata),seqconc(dseq))
		distances <- distances[mcorr]
		names(distances) <- NULL
	}
	else { ## !Refseq
		magicSeq <- order(mcorr)
		magicIndex <- c(unique(rank(mcorr, ties.method="min")), nrow(seqdata)+1)-1
		
		 #SEXP cstringdistanceOO(SEXP Ssequences, SEXP seqdim, SEXP lenS, SEXP paramS, SEXP normS, SEXP magicIndexS, SEXP magicSeqS, SEXP disttypeS)
		
		if(numcpu==1){
			distances <- .Call(TMR_cstringdistanceOO,
				as.integer(dseq),
				as.integer(dim(dseq)),
				as.integer(slength),
				params,
				as.integer(norm),
				as.integer(magicIndex),
				as.integer(magicSeq),
				disttype)
		} else {
			distances <- .Call(TMR_cstringdistanceOOomp,
				as.integer(dseq),
				as.integer(dim(dseq)),
				as.integer(slength),
				params,
				as.integer(norm),
				as.integer(magicIndex),
				as.integer(magicSeq),
				disttype,
				as.integer(numcpu))
		}
		if(method=="NMS"||method=="SVRspell"||method=="NMSMST"){
			distances <- sqrt(distances)
		}
		## Setting some attributes for the dist object
		class(distances) <- "dist"
		attr(distances,"Size") <- length(magicSeq)
		attr(distances,"method") <- method
		attr(distances, "Labels") <- dimnames(seqdata)[[1]]
		attr(distances, "Diag") <- FALSE
		attr(distances, "Upper") <- FALSE
	}

	fin <- proc.time()
	totaltime <- format(round(difftime(as.POSIXct(sum(fin[1:2]), origin="1960-01-01"), as.POSIXct(sum(debut[1:2]), origin="1960-01-01")), 3))
	message(" [>] total time: ", totaltime)

	if (full.matrix && inherits(distances, "dist")) {
		return(TraMineR:::dist2matrix(distances))
	}
	else {
		return(distances)
	}
}
