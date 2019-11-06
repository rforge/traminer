## Table of indicators

seqindic <- function(seqdata, indic=c("visited","trans","ient","cplx"), with.missing=FALSE,
              ipos.args=list(), prec.args=list()) {

	if (!inherits(seqdata,"stslist"))
		stop("data is NOT a sequence object, see seqdef function to create one")

  indic.list <- c("lgth","nonm","dlgth","visited","trans",
    "ntrans","ient","cplx","turb","nturb","all",
    "prec","ipos")
  if (!all(indic %in% indic.list)){
    stop("invalid values in indic: ",paste(indic[!indic %in% indic.list], collapse=", "))
  }

  if (any(indic=="all")) {
    indic.all <- indic.list[-((length(indic.list)-1):length(indic.list))]
    if("ipos" %in% indic) indic.all <- c(indic.all,"ipos")
    if("prec" %in% indic) indic.all <- c(indic.all,"prec")
    indic <- indic.all
  }

  tab <- as.data.frame(rownames(seqdata))
  lab <- "id"

  if("lgth" %in% indic){
  ## Sequence length
    lgth <- suppressMessages(
      seqlength(seqdata, with.missing=TRUE))
    tab <- cbind(tab,lgth)
    lab <- c(lab,"Lgth")
  }
  if("nonm" %in% indic){
  ## Number of non-missing elements
    nonm <- suppressMessages(
      seqlength(seqdata, with.missing=FALSE))
    tab <- cbind(tab,nonm)
    lab <- c(lab,"NonM")
  }
  if("dlgth" %in% indic){
  ## Length of dss
	  dlgth <- suppressMessages(
      seqtransn(seqdata, with.missing=with.missing, norm=FALSE)) + 1
    tab <- cbind(tab,dlgth)
    lab <- c(lab,"Dlgth")
  }
  if("visited" %in% indic){
  ## Number of visited states
    sdist <- suppressMessages(
      seqistatd(seqdata, with.missing=with.missing))
    nvisit <- rowSums(sdist>0)
    tab <- cbind(tab,nvisit)
    lab <- c(lab,"Nvis")
  }
  if("trans" %in% indic){
	## Number of transitions
	  trans <- suppressMessages(
      seqtransn(seqdata, with.missing=with.missing, norm=FALSE))
    tab <- cbind(tab,trans)
    lab <- c(lab,"Ntrn")
  }
  if("ntrans" %in% indic){
	## Proportion of transitions
	  trans <- suppressMessages(
      seqtransn(seqdata, with.missing=with.missing, norm=TRUE))
    tab <- cbind(tab,trans)
    lab <- c(lab,"Ptrn")
  }
  if("ient" %in% indic){
	## Longitudinal Entropy
	  ient <- suppressMessages(
		  seqient(seqdata, with.missing=with.missing, norm=TRUE))
    tab <- cbind(tab,ient)
    lab <- c(lab,"Entr")
  }
  if("cplx" %in% indic){
	## Longitudinal Entropy
	  ici <- suppressMessages(
		  seqici(seqdata, with.missing=with.missing, silent=TRUE))
    tab <- cbind(tab,ici)
    lab <- c(lab,"Cplx")
  }
  if("nturb" %in% indic){
	## Longitudinal Entropy
	  nturb <- suppressMessages(seqST(seqdata, norm=TRUE, with.missing=with.missing, silent=TRUE))
    tab <- cbind(tab,nturb)
    lab <- c(lab,"Turbn")
  }
  if("turb" %in% indic){
	## Longitudinal Entropy
	  turb <- suppressMessages(seqST(seqdata, norm=FALSE, with.missing=with.missing, silent=TRUE))
    tab <- cbind(tab,turb)
    lab <- c(lab,"Turb")
  }

  if("ipos" %in% indic){
  ## index of precarity
    if(!is.null(ipos.args[["seqdata"]])) warning( "[!] seqdata argument in ipos.args will be overwritten!" )
    ipos.args[["seqdata"]] <- seqdata
    if(!is.null(ipos.args[["with.missing"]])) warning( "[!] with.missing argument in ipos.args will be overwritten!" )
    ipos.args[["with.missing"]] <- with.missing
    ipos <- do.call(seqipos, args=ipos.args)
    tab <- cbind(tab,ipos)
    lab <- c(lab,"Ppos")
  }

  if("prec" %in% indic){
  ## index of precarity
    ##!existsFunction("seqprecarity")
    ##  stop("[!] For 'prec', load the TraMineRextras library")
    if(!is.null(prec.args[["seqdata"]]))
      warning( "[!] seqdata argument in prec.args will be overwritten!" )
    prec.args[["seqdata"]] <- seqdata
    if(!is.null(prec.args[["with.missing"]]))
      warning( "[!] with.missing argument in prec.args will be overwritten!" )
    prec.args[["with.missing"]] <- with.missing
    prec <- do.call(seqprecarity, args=prec.args)
    tab <- cbind(tab,prec)
    lab <- c(lab,"Prec")
  }

  names(tab) <- lab
  tab <- tab[,-1, drop=FALSE]
	return(tab)
}
