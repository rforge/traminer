## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2013-12-06
##
## Description:
## Functions that apply to 'vcolmm1' and 'vcolmm2' objects
##
## Overview:
## appendDefArgs:   adds default arguments to a list of arguments
## addEmptyChar:    
## formatMatrix:
## --------------------------------------------------------- #

appendDefArgs <- function(args, default) {
  subs <- setdiff(names(default), names(args))
  if (length(subs) > 0)
    for (i in subs) args[[i]] <- default[[i]] 
  return(args)
}

addEmptyChar <- function(x, nchar, justify = "left") {
  for (i in 1:length(x))
    if (justify == "left") {
      x[i] <- paste(x[i], paste(rep(" ", nchar - nchar(x[i])), collapse = ""), sep = "")
    } else if (justify == "right") {
      x[i] <- paste(paste(rep(" ", nchar - nchar(x[i])), collapse = ""), x[i], sep = "")
    }
  return(x)
}

formatMatrix <- function(x, ...) {
  rowNames <- rownames(x)
  colNames <- colnames(x)
  if (is.null(rowNames)) rowNames <- rep("", nrow(x))
  if (!is.null(colNames)) rowNames <- append("", rowNames)
  rowNames <- addEmptyChar(rowNames, max(nchar(rowNames)))
  rval <- format(x, ...)
  if (!is.null(colNames)) {
    if (max(nchar(c(rval))) < max(nchar(colNames))) {
      for (i in 1:length(rval)) rval[i] <- addEmptyChar(rval[i], max(nchar(colNames)), "right")
    } 
    colNames <- addEmptyChar(colNames, max(nchar(c(rval))), "right")
    rval <- rbind(colNames, rval)
  }
  rval <- cbind(rowNames, rval)
  rval <- apply(rval, 1, paste, collapse = "  ")
  rval <- paste(rval, "\n", sep = "")
  rval <- paste(rval, collapse = "")
  return(rval)
}

tvcolmm_value_space <- function(data, neval = 50L) {

  FUN <- function(x, neval) {
    rval <- sort(unique(x))
    if (is.numeric(x) && length(rval) > neval) {
      rval <- as.double(quantile(rval, seq(0,1,length.out = neval)))
    }
    if (is.factor(x)) rval <- droplevels(rval)
    return(rval)
  }
  return(lapply(as.list(data), FUN, neval = neval))
}
