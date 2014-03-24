## --------------------------------------------------------- #
## Author:          Reto Buergin, rbuergin@gmx.ch
## Date:            2014-03-21
##
## Description:
## Unexported functions imported from other packages
##
## References:
## regr0:           https://r-forge.r-project.org/R/?group_id=555
## lme4:            http://cran.r-project.org/web/packages/lme4/
## party:           http://CRAN.R-project.org/package=party
## strucchange:     http://cran.r-project.org/web/packages/strucchange/
## statmod:         http://cran.r-project.org/web/packages/statmod/
## matrixcalc:      http://cran.r-project.org/web/packages/matrixcalc/
##
## Contents:
## stats;           safe_pchisq
## regr0:           prevgumbel, qrevgumbel, condquant
## lme4:            findbars, nobars, subbars, subnms, slashTerms,
##                  makeInteraction, expandSlash
## party:           logp.supLM
## strucchange:     "sc.maxL2", "sc.meanL2", "sc.beta.sup"
## statmod:         gauss.quad,
## matrixcalc:      vec, u.vectors, E.Matrices, L.matrix
## --------------------------------------------------------- #

## copied from stats:::safe_pchisq
safe_pchisq <- function (q, df, ...) {
    df[df <= 0] <- NA
    pchisq(q = q, df = df, ...)
}

## imports from the regr0 package

nafalse <- function(x) if (is.null(x)) FALSE else ifelse(is.na(x), FALSE, x)

prevgumbel <- function (q, location = 0, scale = 1) {
  if (!nafalse(scale > 0)) 
    stop("\"scale\" must be positive")
  1 - exp(-exp((q - location)/scale))
}


qrevgumbel <- function (p, location = 0, scale = 1) {
  if (!nafalse(scale > 0)) 
    stop("\"scale\" must be positive")
  location + scale * log(-log(p))
}


condquant <- function (x, dist = "normal", sig = 1)  {
  fp <- switch(dist, normal = pnorm, gaussian = pnorm,
               unif = function(x) x, 
               logis = plogis, revgumbel = prevgumbel)
  if (is.null(fp)) 
    stop(paste("!condquant! distribution", dist, "not known"))
  fq <- switch(dist, normal = qnorm, gaussian = qnorm,
               unif = function(x) x, 
               logis = qlogis, revgumbel = qrevgumbel)
  lp <- fp(x)
  lpp <- if (NCOL(lp) == 1) {
    sweep(outer(lp, c(0.5, 0.75, 0.25)), 2, c(0.5, 0.25, 0.75), "+")
  } else {
    lp %*% rbind(c(0.5, 0.75, 0.25), c(0.5, 0.25, 0.75))
  }
  lprand <- lp[, 1] + (lp[, 2] - lp[, 1]) * runif(nrow(lp))
  lr <- cbind(cbind(median = fq(lpp[, 1]), lowq = fq(lpp[, 2]),
                    uppq = fq(lpp[, 3]), random = fq(lprand)) * sig, 
              prob = abs(lp[, 2] - lp[, 1]))
  class(lr) <- "condquant"
  lr
}


## imported from the lme4 package, version 0.999999-0

findbars <- function(term)
### Return the pairs of expressions that separated by vertical bars
{
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
    if (!is.call(term)) stop("term must be of class call")
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
}

## change the case where there is no fixed effect in the model

nobars <- function(term)
### Return the formula omitting the pairs of expressions that are
### separated by vertical bars
{
    if (!('|' %in% all.names(term))) return(term)
    if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
    if (length(term) == 2) {
	nb <- nobars(term[[2]])
	if (is.null(nb)) return(NULL)
	term[[2]] <- nb
	return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

subbars <- function(term)
### Substitute the '+' function for the '|' function
{
    if (is.name(term) || !is.language(term)) return(term)
    if (length(term) == 2) {
	term[[2]] <- subbars(term[[2]])
	return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name('|'))
	term[[1]] <- as.name('+')
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
}

subnms <- function(term, nlist)
### Substitute any names from nlist in term with 1
{
    if (!is.language(term)) return(term)
    if (is.name(term)) {
        if (any(unlist(lapply(nlist, get("=="), term)))) return(1)
        return(term)
    }
    stopifnot(length(term) >= 2)
    for (j in 2:length(term)) term[[j]] <- subnms(term[[j]], nlist)
    term
}

slashTerms <- function(x)
### Return the list of '/'-separated terms in an expression that
### contains slashes
{
    if (!("/" %in% all.names(x))) return(x)
    if (x[[1]] != as.name("/"))
        stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}

makeInteraction <- function(x)
### from a list of length 2 return recursive interaction terms
{
    if (length(x) < 2) return(x)
    trm1 <- makeInteraction(x[[1]])
    trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
    list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
}

expandSlash <- function(bb)
### expand any slashes in the grouping factors returned by findbars
{
    if (!is.list(bb)) return(expandSlash(list(bb)))
    ## I really do mean lapply(unlist(... - unlist returns a
    ## flattened list in this case
    unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
            return(lapply(unlist(makeInteraction(trms)),
                          function(trm) substitute(foo|bar,
                                                   list(foo = x[[2]],
                                                        bar = trm))))
        x
    }))
}


## from package colorspace

"rainbowPalette" <- structure(c("#C3687C", "#AE783C", "#838A0A", "#31964F", "#00998B", "#0090B6", "#847AC4", "#B966AE"))

"divergePalette" <- structure(list(c("#023FA5", "#E2E2E2", "#8E063B"), c("#023FA5", "#BEC1D4", "#D6BCC0", "#8E063B"), c("#023FA5", "#A1A6C8", "#E2E2E2", "#CA9CA4", "#8E063B"), c("#023FA5", "#8C94BF", "#D2D3DC", "#DDD0D2", "#C18692", "#8E063B"), c("#023FA5", "#7D87B9", "#BEC1D4", "#E2E2E2", "#D6BCC0", "#BB7784", "#8E063B"), c("#023FA5", "#727EB5", "#AEB2CD", "#D8D9DE", "#DFD7D9", "#D0AAB1", "#B56B7A", "#8E063B"), c("#023FA5", "#6A76B2", "#A1A6C8", "#CBCDD9", "#E2E2E2", "#DBC9CC", "#CA9CA4", "#B16273", "#8E063B"), c("#023FA5", "#6371AF", "#959CC3", "#BEC1D4", "#DBDCE0", "#E0DBDC", "#D6BCC0", "#C6909A", "#AE5A6D", "#8E063B"), c("#023FA5", "#5D6CAE", "#8C94BF", "#B3B7CF", "#D2D3DC", "#E2E2E2", "#DDD0D2", "#D2B0B6", "#C18692", "#AB5468", "#8E063B")))

## imported from the statmod package, version 1.4.15

gauss.quad <- function(n,kind="legendre",alpha=0,beta=0) {
#	Calculate nodes and weights for Guassian quadrature.
#	Adapted from Netlib routine gaussq.f
#	Gordon Smyth, Walter and Eliza Hall Institute
#	4 Sept 2002. Last modified 6 Aug 2012.

	n <- as.integer(n)
	if(n<0) stop("need non-negative number of nodes")
	if(n==0) return(list(nodes=numeric(0), weights=numeric(0)))
	kind <- match.arg(kind,c("legendre","chebyshev1","chebyshev2","hermite","jacobi","laguerre"))
	i <- 1:n
	i1 <- i[-n] # 1:(n-1)
	switch(kind, legendre={
		muzero <- 2
		a <- rep(0,n)
		b <- i1/sqrt(4*i1^2-1)
	}, chebyshev1={
		muzero <- pi
		a <- rep(0,n)
		b <- rep(0.5,n-1)
		b[1] <- sqrt(0.5)
	}, chebyshev2={
		muzero <- pi/2
		a <- rep(0,n)
		b <- rep(0.5,n-1)
	}, hermite={
		muzero <- sqrt(pi)
		a <- rep(0,n)
		b <- sqrt(i1/2)
	}, jacobi={
		ab <- alpha+beta
		muzero <- exp((ab+1)*log(2) + lgamma(alpha+1) + lgamma(beta+1) - lgamma(ab+2))
		a <- i
		a[1] <- (beta-alpha)/(ab+2)
		i2 <- 2:n
		abi <- ab+2*i2
		a[i2] <- (beta^2-alpha^2)/(abi-2)/abi
		b <- i1
		b[1] <- sqrt(4*(alpha+1)*(beta+1)/(ab+2)^2/(ab+3))
		i2 <- i1[-1] # 2:(n-1)
		abi <- ab+2*i2
		b[i2] <- sqrt(4*i2*(i2+alpha)*(i2+beta)*(i2+ab)/(abi^2-1)/abi^2)
	}, laguerre={
		a <- 2*i-1+alpha
		b <- sqrt(i1*(i1+alpha))
		muzero <- gamma(alpha+1)
	})
	A <- rep(0,n*n)
	A[(n+1)*(i-1)+1] <- a
	A[(n+1)*(i1-1)+2] <- b
	A[(n+1)*i1] <- b
	dim(A) <- c(n,n)
	vd <- eigen(A,symmetric=TRUE)
	w <- rev(as.vector( vd$vectors[1,] ))
	w <- muzero * w^2
	x <- rev( vd$values )
	list(nodes=x,weights=w)
}


## imported from the matrixcalc package, version 1.0-2

vec <- function( x )
{
###
### this function returns a column vector that is a stack of the columns of x
###
### Parameters
### x = a numeric matrix
###
    if ( !is.matrix( x ) ) {
        stop( "argument x is not a matrix" )
    }
    if ( !is.numeric( x ) ) {
        stop( "argument x is not a numeric matrix" )
    }    
    return( t( t( as.vector( x ) ) ) )
}


u.vectors <- function( n )
{
###
### This function constructs an identity matrix I of order
### p = n * ( n +1 ) / 2.  It also builds a lower triangular square matrix of
### order n.  The value of element [i,j] is the column number
### in the identify matrix.  It is the mapping from the coordinates
### to the column vector in the identity matrix.
###
### argument
### n = a positive integer value
###
    if ( n != trunc( n ) )
        stop( "argument n is not an integer" )
    if ( n < 2 )
        stop( "argument n is less than 2" )
    p <- n * ( n + 1 ) / 2
    I <- diag( rep( 1, p ) )
    k <- matrix( 0, nrow=n, ncol=n )
    for ( j in 1:n ) {
        for ( i in j:n ) {
            k[i,j] <- ( j - 1 ) * n + i - 0.5 * j * ( j - 1 )
        }
    }
    return( list( k=k, I=I ) )
}

E.matrices <- function( n )
{
###
### This function creates a list of lists.  The number of components
### in the higher level list is n.  For each component i in the
### higher level list, the number of components in the sub-list
### is n.  Each component j of the sub-list is an n by n matrix
### e_i %*% t ( e_j ).  Each of the arrays is a vector in an
### n by n identity matrix
###
### argument
### n = a positive integer value greater than or equal to 2
###
    if( missing( n ) )
        stop( "argument n is missing" )
    if ( !is.numeric( n ) )
        stop( "argument n is not numeric" )
    if ( n != trunc( n ) )
        stop( "argument n is not an integer" )
    if ( n < 2 )
        stop( "argument n is less than 2" )
    I <- diag( rep( 1, n ) )
    E <- list()
    for ( i in 1:n ) {
        E[[i]] <- list()
        for ( j in 1:n ) {
            E[[i]][[j]] <- I[i,] %o% I[j,]
        }
    }
    return( E )
}

L.matrix <- function( n )
{
###
### This function constructs the elimination matrix as a mapping
### from vec(A) to vech(A)
###
### Arguments
### n = a positive integer value for the order of the matrix
###
    if ( missing( n ) )
        stop( "argument n is missing" )
    if ( !is.numeric( n ) )
        stop( "argument n is not numeric" )
    if ( n != trunc( n ) )
        stop( "argument n is not an integer" )
    if ( n < 2 )
      if ( n == 1) {
        return(matrix(1, 1, 1))
      } else if ( n == 0) {
        return(matrix(, 0, 0))
      } else {
        stop( "argument n is less than 0" )
      }
    u <- u.vectors( n )
    E <- E.matrices( n )
    k <- u$k
    I <- u$I
    p <- n * ( n + 1 ) / 2
    nsq <- n * n
    L <- matrix( 0, nrow=p, ncol=nsq)
    for ( j in 1:n ) {
        for ( i in j:n ) {
            L <- L + I[,k[i,j]] %*% t( vec( E[[i]][[j]] ) )
        }
    }
    return( L )
}
