seqdistOMTrans <- function(seqObj, indel, sm, transmetric="constant", stateweight, norm, with.missing, addcolumn=TRUE, previous=FALSE, ...){
	if(with.missing){
		stop(" [!] with.missing is not supported for this metric!")
	}
	dimnames(sm) <- list(alphabet(seqObj), alphabet(seqObj))
	if(length(indel)==1){
		indel <- rep(indel, length(alphabet(seqObj)))
	}
	names(indel) <- alphabet(seqObj)
	if(transmetric=="prob"){
		tr <- seqtrate(seqObj)
		dimnames(tr) <- list(alphabet(seqObj), alphabet(seqObj))
	}
	seqdata <- as.matrix(seqObj)
	maxcol <- ncol(seqdata)
	newseqdata <-matrix("", nrow=nrow(seqdata), ncol=maxcol)
	sep <- "@@@@TraMineRSep@@@@"

	mypastefunc <- function(i){
		minmax <- function(i){
			return(max(1, min(i, maxcol)))
		}
		ret <- paste( seqdata[, minmax(i)], seqdata[, minmax(i+1)], sep=sep)
		if(previous){
			ret <- paste(seqdata[, minmax(i-1)],ret, sep=sep)
		}
		return(ret)
	}
	
	for (i in 1:(maxcol)) {
		newseqdata[, i] <- mypastefunc(i)
	}
	if(!addcolumn){
		newseqdata <- newseqdata[, -maxcol]
		if(previous){
			newseqdata <- newseqdata[, -1]
		}
	}
	## Setting void states back to NA  (nr will be considered as a distinct state)
	
	alphabet_size <- length(unique(as.character(newseqdata)))
	suppressMessages(newseqdata <- seqdef(newseqdata, cpal=rep("blue", alphabet_size)))
	if(transmetric=="raw"){
		return(seqdistOO(newseqdata, method="OM", sm="CONSTANT", indel=1))
	}
	if(stateweight <0 || stateweight>1){
		stop(" [!] stateweight should be in [0; 1]!")
	}
	transweight <- 1 - stateweight
	if(previous){
		transweight <- transweight/2
	}
	indelrate <- max(indel)/max(sm)
	transweight <- transweight* indelrate
	indel <- indelrate*indel/max(indel)
	
	sm <- sm/max(sm)
	## =========================================
	## Building the new substitution cost matrix
	## =========================================
	
	## Build subsitution matrix and new alphabet
	alphabet <- alphabet(newseqdata)
	alphabet_size <- length(alphabet)
	## Recomputing the subsitution matrix
	indels <- numeric(alphabet_size)
	names(indels) <- alphabet
	newsm <- matrix(0, nrow=alphabet_size, ncol=alphabet_size)
	stateindel <- indel* stateweight
	## indel costs change according to the previous parameters
	if(!previous){
		for(i in 1:alphabet_size){
			statesi <- strsplit(alphabet[i], sep)[[1]]
			indels[i] <- stateindel[statesi[1]]
			if(statesi[1]!= statesi[2]){
				if(transmetric=="constant"){
					 rawtransindel <- 1
				}else if(transmetric=="prob"){
					rawtransindel <- (1-tr[statesi[1], statesi[2]])
				}else if(transmetric=="subcost"){
					rawtransindel <- (sm[statesi[1], statesi[2]])
				}
				indels[i] <- indels[i] +  transweight*rawtransindel
			}
		}
		
	}else{
		for(i in 1:alphabet_size){
			statesi <- strsplit(alphabet[i], sep)[[1]]
			rawtransindel <- 0
			if(transmetric=="constant"){
				 rawtransindel <- (statesi[1]!= statesi[2])+(statesi[2]!= statesi[3])
			}else if(transmetric=="prob"){
				rawtransindel <- (2-tr[statesi[1], statesi[2]]-tr[statesi[2], statesi[3]])
			}else if(transmetric=="subcost"){
				rawtransindel <- (sm[statesi[1], statesi[2]]+sm[statesi[2], statesi[3]])
			}
			indels[i] <- stateindel[statesi[2]] +  transweight*rawtransindel
		}
		
	
	}
	## Just fix the state we are comparing according to 'previous'
	compare_state <- ifelse(previous, 2, 1)
	for (i in 1:(alphabet_size-1)) {
		statesi <- strsplit(alphabet[i], sep)[[1]]
		for (j in (i+1):alphabet_size) {
			statesj <- strsplit(alphabet[j], sep)[[1]]
			cost <- sm[statesi[compare_state], statesj[compare_state]]
			if(transmetric %in% c("constant", "prob", "subcost")){
				cost <- stateweight*cost + (indels[alphabet[i]] + indels[alphabet[j]] -stateindel[statesi[compare_state]] -stateindel[statesj[compare_state]])
			}
			newsm[i, j] <- cost
			newsm[j, i] <- cost
		}
	}
	return(seqdistOO(newseqdata, method="OMloc", sm=newsm, indel=indels, norm=norm, ...))
	
}
