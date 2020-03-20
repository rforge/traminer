seqgranularity <- function(seqdata, tspan=3, method="last"){

  metlist <- c("first","first.valid","last","last.valid","mostfreq")
  if (!method %in% metlist) {
      stop(" [!] method must be one of: ", paste(metlist, collapse = " "),
          call. = FALSE)
  }

  n <- nrow(seqdata)
  missing.char <- attr(seqdata,"missing")
  nr <- attr(seqdata,"nr")
  void <- attr(seqdata, "void")
  lgth <- max(seqlength(seqdata))
  alph <- alphabet(seqdata)

  new.lgth <- ceiling(lgth/tspan)
  new.lgth.f <- floor(lgth/tspan)

  newseq <- seqdata[,1:new.lgth]
  cnames <- names(seqdata)
  newcnames <- cnames[seq(from=1, to=lgth, by=tspan)]

  firstvalid <- function (seq, alphabet, first=TRUE, missing.char=NA) {
  	pos <- which(seq %in% alphabet)
  	if (length(pos)==0) res <- missing.char
  	else {
      if (first) fpos <- min(pos)
      else fpos <- max(pos)
      res <- seq[fpos]
    }
  	return(res)
	}

  prev <- 0
  if (method=="first") prev <- tspan - 1
  for (i in 1:new.lgth.f){
      newseq[,i] <- seqdata[,tspan*i - prev]
  }

  if (method %in% c("first.valid","last.valid")){
    first <- method == "first.valid"
    for (i in 1:new.lgth.f){
      seq <- seqdata[,(tspan*(i-1) + 1):(tspan*i)]
      newseq[,i] <- apply(seq, 1, firstvalid, alphabet=alph, first=first, missing.char=missing.char)
    }
  }
  if (method=="mostfreq") {
      for(i in 1:new.lgth.f) {
          st.freq <- suppressMessages(seqistatd(seqdata[,(tspan*(i-1) + 1):(tspan*i)]))
          newseq[,i] <- apply(st.freq,1,function(x){ifelse(max(x)==0,missing.char,names(which.max(x)))})
      }
  }
  if (new.lgth > new.lgth.f){
      if (method=="first"){
          newseq[,new.lgth] <- seqdata[,tspan*new.lgth.f + 1]
      } else if (method %in% c("first.valid","last.valid")){
          first <- method == "first.valid"
          seq <- seqdata[,(tspan*new.lgth.f+1):lgth]
          newseq[,new.lgth] <- apply(seq, 1, firstvalid, alphabet=alph, first=first, missing.char=missing.char)
      } else if (method=="mostfreq") {
          st.freq <- suppressMessages(seqistatd(seqdata[,(tspan*new.lgth.f+1):lgth]))
          newseq[,new.lgth] <- apply(st.freq,1,function(x){ifelse(max(x)==0,missing.char,names(which.max(x)))})
      } else { # method=="last"
          newseq[,new.lgth] <- seqdata[,lgth]
      }
      newcnames[new.lgth] <- cnames[lgth]
  }
  if (method %in% c("mostfreq","first.valid","last.valid")) {
      newseq <- suppressMessages(seqdef(newseq,
                  alphabet=alph,
                  weights =attr(seqdata,"weights"),
                  labels  =attr(seqdata,"labels"),
                  start   =attr(seqdata,"start"),
                  missing =attr(seqdata,"missing"),
                  nr      =attr(seqdata,"nr"),
                  void    =attr(seqdata,"void"),
                  #xtstep  =attr(seqdata,"xtstep"),
                  cpal    =attr(seqdata,"cpal"),
                  tick.last=attr(seqdata,"tick.last"),
                  Version =attr(seqdata,"Version")
                  ))
  }
  colnames(newseq) <- newcnames
  attr(newseq,"xtstep") <- ceiling(attr(seqdata,"xtstep")/tspan)
	
  return(newseq)
}


#####
## We need to redefine the "[ststlist" method
## will be in TraMineR v2.0-16
## gr

"[.stslist" <- function(x,i,j,drop=FALSE) {
	## Specialized only for column subscript
	## If one column we keep the original data.frame method
	## Otherwise we copy attributes and update "start" value

  ## For negative j, we first build the new subscript set
  if (!missing(j) && j[1]<0) {
    k <- -j
    j <- 1:ncol(x)
    j <- j[! j %in% k]
  }

	if (!missing(j) && length(j)>=1) {
		## Storing the attributes
		x.attributes <- attributes(x)

		## Applying method
	     x <- NextMethod("[")

    if (length(j) == 1) {
      x <- as.data.frame(x)
      class(x) <- c("stslist", class(x))
    }

		## Adapting column names
		x.attributes$names <- x.attributes$names[j]

		## Redefining attributes
		attributes(x) <- x.attributes

	     attr(x,"start") <- x.attributes$start-1+j[1]

		if (!missing(i)) {
			attr(x,"row.names") <- attr(x,"row.names")[i]
			attr(x,"weights") <- attr(x,"weights")[i]
		}

		return(x)
	}

	x <- NextMethod("[")

	if (!missing(i))
		attr(x,"weights") <- attr(x,"weights")[i]

	return(x)
 }
