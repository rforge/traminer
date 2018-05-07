#'@keywords internal
trmatplot.default <- function (d, seed = NULL, 
											rowconstraint = TRUE, morder= 1,
											cspal = NULL, cpal = NULL, title = NULL,
                      xlab =  NULL, ylab = NULL, ylim = NULL, 
											xtlab = NULL, ytlab = NULL,
                      pfilter = NULL,
                      shade.col = "grey80",
                      num = 1,
                      hide.col = NULL,
                      lorder = NULL,
											plot = TRUE,
                      verbose = FALSE, ...){
  
  if ( verbose )
    
    cat ( " [>] check arguments\n" )

	##
	## class  
  if ( ! is.matrix ( d ) ) {
    
    stop ( "[!] d must be a transition probability matrix " ) # is either invalid or missing
    
  }
	
	##
	## matrix
  if ( is.matrix ( d ) ) {
    
    if ( morder == 1 & ! identical ( nrow ( d ) , ncol ( d ) ) ) {
      
      stop ( "[!] if morder is set to equal 1, the transition probability matrix must be square, nrow = ncol" ) 
      
    }
    
    if ( any ( d > 1 ) | any ( d < 0 ) ) {
      
      stop ( "[!] elements of transtion probability matrix must be between zero and one" )
      
    }
    
    if ( rowconstraint == TRUE & sum1 ( d ) == FALSE ){
      
      stop ( "[!] rows of transtion probability matrix must sum up to one, otherwise set argument rowconstraint to FALSE" )
      
    }
    
  }
	
 	##
	## seed
	if ( ! is.null ( seed ) & ! is.numeric ( seed ) ) {
			stop ("[!] seed must be numeric")
	}
	
 	##
	## rowconstraint
	if ( ! is.null ( rowconstraint ) & ! is.logical ( rowconstraint ) ) {
			stop ("[!] rowconstraint must be logical")
	}
	
 	##
	## morder
	if ( ! is.null ( morder ) & ! is.numeric ( morder ) ) {
			stop ("[!] morder must be numeric")
	}
  
	##
	## cspal
  if ( ! is.null (cspal)) {
    
    if (! cspal %in% c ( "dynamic", "harmonic" , "cold", "warm" , "heat" , "terrain" ) ){
      
      stop ( "[!] cspal must be specified as one of: dynamic, harmonic, cold, warm, heat, terrain" )
    }
  }

	##
	## cpal    
  if ( ! is.null ( cpal ) ) {
    
    if ( is.matrix ( d ) ) {
      
      if ( length ( cpal ) != dim ( d ) [ 1 ] ) {
        
        stop("[!] number of colors in 'cpal' must equal number of states")
      }
    }
    if ( ! is.null ( cspal ) ){
      
      stop ( "[!] only one of cpal or cspal can be specified as non-null" )
      
    }
  }
         
  
 	##
	## pfilter 
  
  if ( ! is.null ( pfilter ) ) {
    
    if  ( ! (pfilter %in% c ( "smax", "smin" ,"tmax", "tmin" ) ) ) {
      
      stop ( "[!] pfilter must be specified as one of: smax, smin, tmax or tmin" )
      
    }
    
    if ( ! is.numeric ( num ) ){
      
      stop ( "[!] num must be numeric")
      
    }
    
    
    if ( is.matrix ( d ) ) {
      
      if ( ! is.element ( num, c( 0 : (dim ( d ) [ 1 ] * dim ( d ) [ 2 ] ) - 1 ) ) ) {
        
        stop("[!] num must be an element of the set (0: matrix size-1)")
        
      }
      
    }
    
  }
  
  
  
  
 	set.seed(seed)
  
 	M <- dim ( d ) [ 2 ]
    
 	l <- morder 
    
	w <- d
 
    
  
  p <- paths ( M , l )
  
  w <- matrix ( t ( w ) , ncol = 1)
  
  w <- as.vector ( w )
  
  s <- suppressMessages ( seqdef ( p , weights = w) )
  
  st <- apply ( p , 1 , function ( x ) paste ( x , collapse = "-" ) )
  
  if ( is.null ( cpal ) ) {
    
    if ( is.null ( cspal ) ) {
      
      cpl <- rainbow_hcl ( M , c = 80 , l = 65 , start = 0 , end = 360 * ( M - 1 ) / M ) # cite: colorspace:::rainbow_hcl
    }
    
    else if ( ! is.null ( cspal ) & cspal == "dynamic"){
      
      cpl <- rainbow_hcl ( M , c = 80 , l = 60 , start = 120, end = 390 ) # cite: colorspace:::rainbow_hcl 
    }
    
    else if ( ! is.null ( cspal ) & cspal == "harmonic"){
      
      cpl <- rainbow_hcl ( M , c = 80 , l = 60 , start = 60 , end = 240 ) # cite: colorspace:::rainbow_hcl    
    }
    
    else if ( ! is.null ( cspal ) & cspal == "cold"){
      
      cpl <- rainbow_hcl ( M , c = 80 , l = 60 , start = 270 , end = 150 ) # cite: colorspace:::rainbow_hcl
    }
    
    else if ( ! is.null ( cspal ) & cspal == "warm"){
      
      cpl <-rainbow_hcl ( M , c = 80 , l = 60 , start = 90 , end = -30 ) # cite: colorspace:::rainbow_hcl
    }
    
    else if ( ! is.null ( cspal ) & cspal == "heat"){
      
      cpl <- heat_hcl ( M , h = c ( 0, 90 ), c. = c ( 100, 100 ), l = 50, power = c ( 1/M, 1 ) )  # cite: colorspace:::heat_hcl
    }
  
    else if ( ! is.null ( cspal ) & cspal == "terrain"){
    
      cpl <- terrain_hcl ( M, h = c ( 130,  0 ), c = c ( 100, 30 ), l = 50, power = c ( 1/M, 1 ) ) # cite: colorspace:::terrain_hcl
    }
  }
  
  else {
    cpl <- cpal
  }
  
  ch <- rep( c ( cpl ) , M ^ ( l - 1 ) , each = M ) # color from first time period
  
  #ch <- rep( c(cpl), M ^ (l)) # color from last time period
  
  predat <- data.frame ( w = w , ch = ch , s = st )
  
  predat <- predat [ order ( predat $ s ) , ] # order by seq name, i.e. state
  
  if ( is.null ( pfilter ) ) {
    
    predat <- predat
    
  }
  
  else if ( ! is.null ( pfilter ) & pfilter == "smin" ) {
    
    predat <- ismin ( M , l , dt = predat)
    
    predat <- smin ( M , l , dt = predat , shade.col )
    
  }
  
  else if ( ! is.null ( pfilter ) & pfilter == "smax" ) {
    
    predat <- ismax ( M , l , dt = predat)
    
    predat <- smax ( M , l , dt = predat , shade.col )
    
  }
  
  else if ( ! is.null ( pfilter ) & pfilter == "tmax" )  {
    
    predat <- tmax ( M , l , dt = predat , shade.col , num = num )
    
  }
  
  else if ( ! is.null ( pfilter ) & pfilter == "tmin" ) {
    
    predat <- tmin ( M , l , dt = predat , shade.col , num = num )
    
  }
  
  dat <- data.frame ( w = as.vector ( predat $ w ) ,
                      ch = as.vector ( predat $ ch ) , 
                      s = as.vector ( predat $ s ) )
  
  dat <- dat [ order ( dat $ s ) , ] # order by seq name
  
  dat <- dat [ order ( dat $ w , decreasing = TRUE ) , ]
  
  # TITLE
  
  if ( is.null ( title ) ) {
    
    ttl <- "Probability Transition Matrix"
    
  }
  
  else if ( ! is.null ( title ) ) {
    
    ttl <- title
    
  }
  
  
  # title
  
  if ( is.null ( xlab ) ) {
    
    xlb <- "Time"
    
  }
  
  else if ( ! is.null ( xlab ) ) {
    
    xlb <- xlab
    
  }
  
  # ylab
  
  if ( is.null(ylab) ){
    
    ylb <- "States"
    
  }
  
  else if ( ! is.null ( ylab ) ) {
    
    ylb <- ylab
    
  }
  
  # ylim    
  
  if ( is.null ( ylim ) ){
    
    ylm <- c ( 0.5 , ( M + 0.5 ) )
    
  }
  
  else if ( ! is.null ( ylim ) ){
    
    ylm <- ylim
    
  }
  
  
  # xtlab
  
  if ( is.null ( xtlab ) ) {
    #xt <- c(0:(M-l))
    #xt <- c(c((-l):0)) # number backwards
    
    xt <- c ( c ( 1 : l ) ) # number forwards
    
    xt <- paste( "t +" , xt )
    
    xt <- c ( "t", xt )
    
  }
  
  else if ( ! is.null ( xtlab ) ) {
    
    xt <- xtlab
    
  }
  
  # alphabet: labeling the y-axis ticks with the visible states (or why not even the hidden states?)
  # yt <-d @ y @ dictionary
  
  # hide.col
  
  if ( ! is.null ( pfilter ) & ! is.null ( hide.col ) ){
    
    hd.col <- hide.col
    
  }
  
  else {
    
    hide.col <- shade.col
    
  }
  
  
  # foreground / background
  
  if ( ! is.null ( pfilter ) ) {
    
    if ( pfilter %in% c ( "smax", "tmax" ) ){
      
      lordr <- "foreground"
      
    }
    
    else if ( pfilter %in% c ( "smin" , "tmin" ) ) {
      
      lordr<-"background"
      
    }
	}
    
  else if ( is.null ( pfilter ) ) {
    
    lordr <- NULL
    
  }

	#	seqpcplot	
	a <- seqpcplot ( seqdata = s, title = ttl, ylab = ylb, xlab = xlb, hide.col = hide.col, lorder = lordr,
         order.align="time", ylim = ylm, cpal= dat$ch, xtlab = xt, verbose = verbose, plot = FALSE, ...) 

	# ytlab

	if ( ! is.null ( ytlab ) ){
		a$ylevs <- ytlab
	}
	else {
		a$ylevs
	}

	#	list arguments inhert in trmatplot
	b <- structure ( list ( rowconstraint = rowconstraint, 
													cspal = cspal, 
													ytlab = ytlab,
													pfilter = pfilter, 
													shade.col = shade.col, 
													num = num ),
									class = "trmatplot")

	#	aggregate list of arguments
	rval <- structure ( list ( plot = a, trmatplot = b, seed = seed, verbose = verbose), class = "trmatplot")

	#	plot

	if ( plot == TRUE ) {

		plot ( a )

	}
	
	else {

		return ( rval )

	}
	
	## DATA
  invisible ( rval )
}
