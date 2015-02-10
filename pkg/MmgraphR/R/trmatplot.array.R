#'@keywords internal
trmatplot.array <- function (d, rowconstraint = TRUE, seed = NULL,
													cspal = NULL, cpal = NULL, title = NULL,
                	       	xlab =  NULL, ylab = NULL, ylim = NULL, xtlab = NULL,
                      		pfilter = NULL,
                       		shade.col = "grey80",
                       		num = 1,
                       		hide.col = NULL,
                       		lorder = NULL,
                       		verbose = FALSE, ...){

	##CHECK
  # M by M matrix
#	if ( dim ( d ) [ 1 ] == dim ( d ) [ 2 ] )
  # M by M^l matrix

#	if ( as.integer ( log ( dim(array, base= 2))==(log(8,base=2))


  d <- matrix ( d, nrow = ( dim ( d )[ 1 ] ), ncol = ( dim ( d )[ 2 ] ))
 
  trmatplot.default ( d = d, rowconstraint = rowconstraint, seed = seed,
											cspal = cspal, cpal = cpal, title = title, 
											xlab = xlab, ylab = ylab, ylim = ylim, xtlab = xtlab,
                      pfilter = pfilter,
                      shade.col = shade.col,
                      num = num,
                      hide.col = hide.col,
                      lorder = lorder,
                      verbose = verbose, ... )

}
