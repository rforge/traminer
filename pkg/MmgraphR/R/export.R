##' -------------------------------------------------------- #
##' Author:          Pauline (Poulcheria) Adamopoulou, padamopo@gmail.com
##' Date:            2014-12-03
##'
##' Description:
##' Additional functions that may be exported
##'
##' Contents:
##'
##' march2traminer  :    converts march.Dcmm class objects to
##' traminer2march  :    convert TraMineR 
##' -------------------------------------------------------- #

##
# march2traminer
##
march2traminer <- function ( d, weights = TRUE, id = TRUE, cnames = TRUE, cpal = NULL, verbose = TRUE,... ) {
  
  if ( weights == TRUE ) {
    m2t <- seqdef ( d @ yRaw, weights = d @ weights)
		if (verbose)
	    cat ( " [>] data transfered\n" )  
		if (verbose)
	    cat ( " [>] weights extracted\n" )    
  }
  
  else {
    m2t <- seqdef ( d @ yRaw )
		if (verbose)
	    cat ( " [>] data transfered\n" ) 
		if (verbose) 
  	  cat ( " [>] weights extracted\n" )    
  }

  if ( !is.null ( cpal ) ) {

    if ( cpal == "specified" ) {
      attr ( m2t, "cpal" ) <- names ( d @ dictionary )
    }

    else {
      attr ( m2t, "cpal" ) <- cpal
    }

  }

  if ( id == TRUE ) {
    attr ( m2t, "row.names" ) <- rownames ( d @ yRaw )
  }

  if ( cnames == TRUE ) {
    attr ( m2t, "names" ) <- colnames ( d @ yRaw )
  }

  m2t
}

#for checks (verbose==TRUE)
#if (length(cpal) != nbstates) 
#stop("\n [!] number of colors in 'cpal' must equal length of alphabet",  call. = FALSE)

##
# traminer2march
##
traminer2march <- function ( d, weights = TRUE, id = TRUE, cnames = TRUE, cpal = TRUE, verbose = TRUE,... ) {

  if ( weights == TRUE & !is.null ( attr ( d, "weights") ) ) {
    ww <- attr ( d, "weights" )
    t2m <- march.dataset.loadFromDataFrame ( d, weights = ww )
		if (verbose)
			cat ( " [>] data transfered\n" )  
		if (verbose)
			cat ( " [>] weights extracted\n" )      
  }

  else if ( weights == FALSE | is.null( attr ( d, "weights") ) ){
    t2m <- march.dataset.loadFromDataFrame ( d )
		if (verbose)
    	cat ( " [>] data transfered\n" )  
		if (verbose)
    	cat ( " [>] no weights specified, none extracted\n" ) 
  }

  if ( id == TRUE ){
    rn <- attr ( d, "row.names" )
    rownames (t2m @ yRaw ) <- rn
    names ( t2m @ y ) <- rn
    names ( t2m @ weights ) <- rn
		if (verbose)
	    cat ( " [>] id extracted\n" ) 
  }

  if ( cnames == TRUE ){
    cn <- attr( d, "names" )
    colnames ( t2m @ yRaw ) <- cn
		if (verbose)
  	  cat ( " [>] column names extracted\n" ) 
  }

  if ( cpal == TRUE ){
    cp <- attr ( d, "cpal" )
    names ( t2m @ dictionary ) <- cp
		if (verbose)
			cat ( " [>] color palette extracted\n" ) 
  }

  else {
		if (verbose)
    	cat ( " [>] no color palette extracted\n" ) 
  }

  t2m
}
