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
march2traminer <- function ( data, weights = TRUE, id = TRUE, cnames = TRUE, cpal = NULL, ... ) {
  
  if ( weights == TRUE ) {
    m2t <- seqdef ( data @ yRaw, weights = data @ weights)
    cat ( " [>] data transfered\n" )  
    cat ( " [>] weights extracted\n" )    
  }
  
  else {
    m2t <- seqdef ( data @ yRaw )
    cat ( " [>] data transfered\n" )  
    cat ( " [>] weights extracted\n" )    
  }

  if ( !is.null ( cpal ) ) {

    if ( cpal == "specified" ) {
      attr ( m2t, "cpal" ) <- names ( data @ dictionary )
    }

    else {
      attr ( m2t, "cpal" ) <- cpal
    }

  }

  if ( id == TRUE ) {
    attr ( m2t, "row.names" ) <- rownames ( data @ yRaw )
  }

  if ( cnames == TRUE ) {
    attr ( m2t, "names" ) <- colnames ( data @ yRaw )
  }

  m2t
}

#for checks (verbose==TRUE)
#if (length(cpal) != nbstates) 
#stop("\n [!] number of colors in 'cpal' must equal length of alphabet",  call. = FALSE)

##
# traminer2march
##
traminer2march <- function ( data, weights = TRUE, id = TRUE, cnames = TRUE, cpal = TRUE, ... ) {

  if ( weights == TRUE ) {
    ww <- attr ( data, "weights" )
    t2m <- march.dataset.loadFromDataFrame ( data, weights = ww )
    cat ( " [>] data transfered\n" )  
    cat ( " [>] weights extracted\n" )      
  }

  else {
    t2m <- march.dataset.loadFromDataFrame ( data )
    cat ( " [>] data transfered\n" )  
    cat ( " [>] no weights specified, none extracted\n" ) 
  }

  if ( id == TRUE ){
    rn <- attr ( data, "row.names" )
    rownames (t2m @ yRaw ) <- rn
    names ( t2m @ y ) <- rn
    names ( t2m @ weights ) <- rn
    cat ( " [>] id extracted\n" ) 
  }

  if ( cnames == TRUE ){
    cn <- attr( data, "names" )
    colnames ( t2m @ yRaw ) <- cn
    cat ( " [>] column names extracted\n" ) 
  }

  if ( cpal == TRUE ){
    cp <- attr ( data, "cpal" )
    names ( t2m @ dictionary ) <- cp
    cat ( " [>] color palette extracted\n" ) 
  }

  else {
    cat ( " [>] no color palette extracted\n" ) 
  }

  t2m
}
