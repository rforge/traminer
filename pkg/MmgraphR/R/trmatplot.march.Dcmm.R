# renamed from trmatplot.march.Dcmm b/c implies S3 method (impossible since march isn't on CRAN)
march.Dcmm.trmatplot <- function(d, seed = NULL, type = "hidden", hstate = 1, 
                       	cspal = NULL, cpal = NULL, main = NULL,
                      	xlab =  NULL, ylab = NULL, ylim = NULL, 
												xtlab = NULL, ytlab = NULL,
												pfilter = NULL,                       	
												shade.col = "grey80",
                       	num = NULL,
                       	hide.col = NULL,
                       	lorder = NULL,
												plot = TRUE,
                       	verbose = FALSE, ...) {

	###
  ### Check Arguments
	###
  if (verbose)
			cat("[>] initial argument check\n")
  ##
	## class
	if (!class(d)[1] == "march.Dcmm") {
			stop("[!] input must be an object of class march.Dcmm")
	}

	##
	## seed
	if (!is.null(seed) & !is.numeric(seed)) {
			stop("[!] seed must be numeric")
	}

  ##
  ## type
  if (!(type %in% c("hidden", "visible"))) {
      stop("[!] matrix type must be specified as either hidden or visible")
  }
	if (type == "visible" & d@orderVC == 0) {
      stop ("[!] visible matrix cannot be plotted when the visible order is zero")
  }

 	##
  ## hstate
  if (!is.element(hstate, c(1:(d@M)))) {
      stop("[!] hstate must be an element of the set of hidden states")
  }

  ##  
  ## cpal
  if (!is.null(cpal)) {
    if (type == "hidden" & length(cpal) != d@M) {
			stop("[!] number of colors in 'cpal' must equal number of hidden states")
    } else if (type == "visible" & length(cpal) != d@y@K) {
			stop("[!] number of colors in 'cpal' must equal number of visible states")
    }
	}
	
	##
  ## cspal  
  if (!is.null(cspal)) {
		if (length(cspal) > 1 | !cspal %in% c("dynamic", "harmonic", "cold", "warm", "heat", "terrain")) {
			stop("[!] cspal is an argument of length one and must be specified as one of: dynamic, harmonic, cold, warm, heat, terrain")
    }
    if (!is.null(cpal)) {
      stop("[!] only one of cpal or cspal can be specified as non-null")
		}
	}
      
	##
	## pfilter 
	if (is.null(pfilter)) {
		if (!is.null(num)) {
			stop("[!] num should be left as null. num only needs to be specified when pfilter is specified as either tmax or tmin")	
		}
	} else { # !is.null(pfilter)
		if (length(pfilter) > 1) {
			stop("[!] pfilter is an argument of length one. If non-null pfilter must be specified one of: smax, smin, tmax or tmin")
		} 
		
		if (!pfilter %in% c("smax", "smin", "tmax", "tmin")) {
			stop("[!] pfilter must be specified as one of: smax, smin, tmax or tmin")
		}
		# num  
		if (!pfilter %in% c("tmax", "tmin")) {
			if (!is.null(num)) {
				stop("[!] num should be left as null. num only needs to be specified when pfilter is specified as either tmax or tmin")		
			}	
		} else if (pfilter %in% c("tmax", "tmin")) {
			if (is.null(num)) {
				stop(paste("[!] num must be specified as a whole number, in this case between 1 and ", length(d) - 1, " when pfilter is specified as either tmax or tmin",sep = ""))
			} else {
				if (length(num) > 1 | !is.numeric(num)) {
  				stop(paste("[!] num is a numeric argument of length one. num must specified as a whole number,in this case between 1 and ", length(d) - 1, " when pfilter is specified as either tmax or tmin", sep = ""))
  			} 
				if (type == "visible") {
					if (is.null(num) | !is.element(num, c(0:(d@y@K ^ (d@orderVC + 1))))) {
    				stop(paste("[!] num must be specified as a whole number, in this case, between 0 and ", d@y@K ^ (d@orderVC + 1), " when pfilter is specified as either tmax or tmin", sep = ""))
					}
				} else if (type == "hidden") {
					if (is.null(num) | !is.element (num, c(0:(d@M ^ (d@orderHC + 1))))) {
    				stop(paste("[!] num must be specified as a whole number, in this case, between 0 and ", d@M ^ (d@orderHC + 1), " when pfilter is specified as either tmax or tmin", sep = ""))
					}
				}
			}
		}
	} 



	###
  ### Create matrix
	###
  set.seed(seed)

  if (type == "hidden") {
    M <- d@M
    l <- d@orderHC
    w <- march.dcmm.h.compactA(d)
  } else if (type == "visible") {
 		M <-  d@y@K # K
    l <- d@orderVC # f
    w <- d@RB[hstate, , ]
  }
  

	###
  ### BEGIN COMPUTATION
	###

  p <- paths(M, l)
  w <- matrix(t(w), ncol = 1)
  w <- as.vector(w)
  s <- suppressMessages(seqdef(p, weights = w))
  st <- apply(p, 1, function(x) paste(x, collapse = "-"))
  
	if (is.null(cpal)) {  
		if (is.null(cspal)) {
			# if (type == "hidden") {
      cpl <- rainbow_hcl(M, c = 80, l = 65, start = 0, end = 360 * (M - 1) / M) # cite: colorspace:::rainbow_hcl
			# } else	if (type == "visible") {
      #		cpl <- heat_hcl(M, h = c(100, -100), l = c(75, 40), c = c(40, 80), power = 8) # ccolorspace:::heat_hcl
      # }
		} else {
			if (cspal == "dynamic") {
				cpl <- rainbow_hcl (M, c = 80, l = 60, start = 120, end = 390) # cite: colorspace:::rainbow_hcl
			} else if (cspal == "harmonic") {
      	cpl <- rainbow_hcl (M, c = 80, l = 60, start = 60, end = 240) # cite: colorspace:::rainbow_hcl    
    	} else if (cspal == "cold") {
	      cpl <- rainbow_hcl(M, c = 80, l = 60, start = 270, end = 150) # cite: colorspace:::rainbow_hcl
    	} else if (cspal == "warm") {
	      cpl <-rainbow_hcl(M, c = 80, l = 60, start = 90, end = -30) # cite: colorspace:::rainbow_hcl
    	} else if (cspal == "heat") {
	      cpl <- heat_hcl(M, h = c(0, 90), c. = c(100, 100), l = 50, power = c(1 / M, 1))  # cite: colorspace:::heat_hcl
    	}	else if (cspal == "terrain") {    
		    cpl <- terrain_hcl(M, h = c(130, 0), c = c(100, 30), l = 50, power = c(1 / M, 1)) # cite: colorspace:::terrain_hcl
    	}
  	}
	} else {
    cpl <- cpal
  }
  
  ch <- rep(c(cpl), M ^ (l - 1), each = M) # color from first time period
  #ch <- rep( c(cpl), M ^ (l)) # color from last time period
  predat <- data.frame(w = w, ch = ch, s = st)
  predat <- predat[order(predat$s), ] # order by seq name
  
  if (is.null(pfilter)) {
		predat <- predat
  } else if (!is.null(pfilter)) {
		if (pfilter == "smin") {
    	predat <- is.min(M, l, dt = predat)
    	predat <- smin(M, l, dt = predat, shade.col)
    } else if (pfilter == "smax") {
  	  predat <- is.max(M, l, dt = predat)
  	  predat <- smax(M, l, dt = predat, shade.col)
    } else if (pfilter == "tmax") {    
    	predat <- tmax(M, l, dt = predat, shade.col, num = num)
    } else if (!is.null(pfilter) & pfilter == "tmin") {
      predat <- tmin(M, l, dt = predat, shade.col, num = num)
    }
	}
  
  dat <- data.frame(w = as.vector(predat$w),
                    ch = as.vector(predat$ch), 
                    s = as.vector(predat$s))
  dat <- dat[order(dat$s), ] # order by seq name
  dat <- dat[order(dat$w, decreasing = TRUE), ]
  
  # TITLE
  if (is.null(main)) {
		if (type == "hidden") {      
      ttl <- "Hidden Transition Matrix"
    } else if (type == "visible") {      
      ttl <- paste("Visible Transition Matrix for Hidden State", hstate)
    } else {		
			ttl <- "Probability Transition Matrix"
		}
  } else { 
    ttl <- main
  }
    
  # XLAB
  if (is.null(xlab)) {
      xlb <- "Time"
  } else {
    xlb <- xlab
  }
  
  # YLAB
  if (is.null(ylab)) {
  	if ( type == "hidden" ) {
  	  ylb <- "Hidden States"      
    } else if ( type == "visible" ) {
      ylb <- "Visible States"
    } else {
	  	ylb <- "States"
    }
  } else {
    ylb <- ylab
  }
  
  #YLIM    
  if (is.null(ylim)) {
		ylm <- c(0.5, (M + 0.5))
  } else if (!is.null(ylim)) {
	  ylm <- ylim
  }
  
  # XtLAB
  if (is.null(xtlab)) {
    #xt <- c(0:(M-l))
    #xt <- c(c((-l):0)) # number backwards
    xt <- c ( c ( 1 : l ) ) # number forwards
    xt <- paste( "t +" , xt )
    xt <- c ( "t", xt )
  } else {    
    xt <- xtlab
  }
  
  # hide.col
  if (!is.null(pfilter) & !is.null(hide.col)) {
  	hd.col <- hide.col
  } else {   
    hide.col <- shade.col
  }
  
	# foreground / background 
  if (!is.null(pfilter)) {    
    if (pfilter %in% c("smax", "tmax")) {      
      lordr <- "foreground"
     } else if (pfilter %in% c ("smin", "tmin")) {
      lordr<-"background"
    }
  } else {
	  lordr <- lorder
	}

	#	seqpcplot	
	a <- seqpcplot(seqdata = s, main = ttl, ylab = ylb, xlab = xlb, 
								hide.col = hide.col, lorder = lordr, order.align = "time",
								ylim = ylm, cpal = dat$ch, xtlab = xt, verbose = verbose,
								plot = FALSE, ...) 

	# ytlab
  # alphabet: labeling the y-axis ticks with the visible states (or why not even the hidden states?)
  # yt <-d @ y @ dictionary
	if (!is.null(ytlab)) {
		a$ylevs <- c(ytlab, "")
	}
	else {
		a$ylevs
	}

	#	list arguments inhert in trmatplot
	b <- structure(list(type = type, 
											hstate = hstate, 
											cspal = cspal, 
											ytlab = ytlab,
											pfilter = pfilter, 
											shade.col = shade.col, 
											num = num), class = "trmatplot")

	#	aggregate list of arguments
	rval <- structure(list(plot = a, 
												 trmatplot = b,
												 seed = seed,
												 verbose = verbose), class = "trmatplot")

	#	plot
	if (plot == TRUE) {
		plot (a)
	} else {
		return(rval)
	}
	
	## DATA
  invisible(rval)
}

