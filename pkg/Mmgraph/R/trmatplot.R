#'@title Transition Matrix Plot
#'@aliases trmatplot
#'@description A coordinate plot which maps each row in the probability transition matrix as a line, where each line is weighted by probability. Users can apply simple filters to empasize the most (or least) probable state sequences overall, or by initial state.
#'@author Poulcheria Adamopoulou
#'@param d Object to be plotted. A transition matrix of probabilities. 
#'@param seed A single value, interpreted as an integer, or \code{NULL} (default). See 'Details'.
#'@param cpal Color palette vector when coloring probability sequences. The \command{\link[colorspace]{rainbow_hcl}} command is used to generate  a color palette if none is specified.
#'@param title Title for the graphic. Default is \emph{'Probability Transition Matrix'}.
#'@param xlab Label for the x axis. Default is \emph{'Time'}.
#'@param ylab Label for the y axis. Default is \emph{'States'}.
#'@param ylim Numeric vector of length 2 giving the y coordinates range.
#'@param xtlab Label for the x axis ticks. Default is time (t, t+1,...).
#'@param pfilter Probability filter. Can be specified as one of \code{"tmax", "tmin", "smax", "smin"}. See 'Details'.
#'@param shade.col The color for sequences shaded out using the \code{pfilter} argument. Default is \code{"grey80"}. See 'Details'.
#'@param num Numeric. The number of sequences to be highlighted. Only applicable when using \code{pfilter="tmax"} or \code{pfilter="tmin"}. Default is 1.
#'@param hide.col The color for sequences shaded out using the \code{filter} argument. Default is \code{"grey80"}. See 'Details'.
#'@param lorder Line order. Either \code{"background"} or \code{"foreground"}. When \code{pfilter} is used \code{lorder} is set by default.
#'@param verbose Logical. Reports extra information on progress. Default is \code{FALSE}.
#'@param ... arguments, such as graphical parameters, to be passed on. See \command{\link[graphics]{par}} and \command{\link[TraMineR]{seqpcplot}}.
#'@details 
#'\describe{ 
#'  Setting a \code{seed} allows the graphic to be replicated.
#'} 
#'\describe{
#'The \code{pfilter} argument serves to highlight probability sequences that are either most probable while shading out those that are less probable in \code{shade.col} and vice-versa. 
#'The four options for \code{pfilter} are described below, see also 'Examples'.
#'\item{\code{"smax"}}{For each initial state the most probable next state is highlighed.}
#'\item{\code{"smin"}}{For each initial state the least probable next state is highlighed.} 
#'\item{\code{"tmax"}}{The sequence of states with the highest probability overall is highlighed. To highlight the \emph{n} most probable sequences of states, set \code{num = n}.}
#'\item{\code{"tmin"}}{The sequence of states with the lowest probability overall is highlighed. To highlight the \emph{n} least probable sequences of states, set \code{num = n}.}
#'}
#'\describe{
#'  The \code{filter} and \code{hide.col} are arguments inherent in and may be passed on to \command{\link[TraMineR]{seqpcplot}}. 
#'  The \code{filter} argument serves to specify filters to gray less interesting patterns. 
#'  The filtered-out patterns are displayed in the \code{hide.col} color. 
#'  The filter argument expects a list with at least elements \code{type} and \code{value}. 
#'  Most relevant within the context of probabilities is \code{type = "sequence"}, which highlights the specific pattern, for example assign 
#'  \code{filter = list(type = "sequence", value = "(1)-(4)")}. 
#'  See 'Examples'.
#'}
#'@return \code{trmatplot} returns an object of class \command{\link[TraMineR]{seqpcplot}}
#'@references  B\"{u}rgin, R. and G. Ritschard (2014), "A decorated parallel coordinate plot for categorical longitudinal data", \emph{The American Statistician}. Vol. 68(2), pp. 98-103. 
#'@seealso \command{\link[TraMineR]{seqpcplot}}
#'@example examples/trmatplot_example.R
#'@export
#'@importFrom colorspace rainbow_hcl
#'@importFrom TraMineR seqdef seqpcplot
trmatplot <- function (d, seed = NULL, cpal = NULL , title = NULL,
                       xlab =  NULL, ylab = NULL, ylim = NULL, xtlab = NULL,
                       pfilter = NULL,
                       shade.col = "grey80",
                       num = 1,
                       hide.col = NULL,
                       lorder = NULL,
                       verbose = FALSE, ...){
  
  if ( verbose )
    
    cat ( " [>] check arguments\n" )
  
  if ( ! is.matrix ( d ) ) {
    
    stop ( "[!] d must be a transition probability matrix " ) # is either invalid or missing
    
  }
  if ( is.matrix ( d ) ) {
    
    if ( ! identical ( nrow ( d ) , ncol ( d ) ) ) {
      
      stop ( "[!] the transition probability matrix must be square, nrow = ncol" ) 
      
    }
    
    if ( any ( d > 1 ) | any ( d < 0 ) ) {
      
      stop ( "[!] elements of transtion probability matrix must be between zero and one" )
      
    }
    
    if ( sum1 ( d ) == FALSE ){
      
      stop ( "[!] rows of transtion probability matrix must sum up to one" )
      
    }
    
  }
  
  
  if ( ! is.null ( cpal ) ) {
    
    if ( is.matrix ( d ) ) {
      
      if ( length ( cpal ) != dim ( d ) [ 1 ] ) {
        
        stop("[!] number of colors in 'cpal' must equal number of states")
      }
    }
  }         
  
  
  
  if ( ! is.null ( pfilter ) ) {
    
    if  ( ! (pfilter %in% c ( "smax", "smin" ,"tmax", "tmin" ) ) ) {
      
      stop ( "[!] pfilter must be specified as either smax, smin, tmax or tmin" )
      
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
  
  
  if ( is.matrix ( d ) ) {
    
    M <- dim ( d ) [ 1 ]
    
    l <- 1 
    
    w <- d
    
  }
  
  
  p <- paths ( M , l )
  
  w <- matrix ( t ( w ) , ncol = 1)
  
  w <- as.vector ( w )
  
  s <- seqdef ( p , weights = w)
  
  st <- apply ( p , 1 , function ( x ) paste ( x , collapse = "-" ) )
  
  if ( is.null ( cpal ) ) {
    
    cpl <- rainbow_hcl ( M , c = 80 , l = 65 , start = 0 , end = 360 * ( M - 1 ) / M ) # cite: rainbow_hcl {colorspace}
    
  }
  
  else {
    
    cpl <- cpal
    
  }
  
  ch <- rep( c ( cpl ) , M ^ ( l - 1 ) , each = M ) # color from first time period
  
  #ch <- rep( c(cpl), M ^ (l)) # color from last time period
  
  predat <- data.frame ( w = w , ch = ch , s = st )
  
  predat <- predat [ order ( predat $ s ) , ] # order by seq name
  
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
  
  
  # XLAB
  
  if ( is.null ( xlab ) ) {
    
    xlb <- "Time"
    
  }
  
  else if ( ! is.null ( xlab ) ) {
    
    xlb <- xlab
    
  }
  
  # YLAB
  
  if ( is.null(ylab) ){
    
    ylb <- "States"
    
  }
  
  else if ( ! is.null ( ylab ) ) {
    
    ylb <- ylab
    
  }
  
  #YLIM    
  
  if ( is.null ( ylim ) ){
    
    ylm <- c ( 0.5 , ( M + 0.5 ) )
    
  }
  
  else if ( ! is.null ( ylim ) ){
    
    ylm <- ylim
    
  }
  
  
  # XtLAB
  
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
  
  if ( ! is.null ( filter ) & ! is.null ( hide.col ) ){
    
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
    
    lordr <- lorder
    
  }
  
  
  seqpcplot(seqdata = s, title = ttl, ylab = ylb, xlab = xlb, hide.col = hide.col, lorder = lordr,
            order.align="time", ylim = ylm, cpal= dat$ch, xtlab = xt, verbose = verbose, ...) #
  
}
