#'@keywords internal
trmatplot.array <- function (d, seed = NULL, 
													rowconstraint = TRUE, morder = 1,
													cspal = NULL, cpal = NULL, main = NULL,
                	       	xlab =  NULL, ylab = NULL, ylim = NULL, 
													xtlab = NULL, ytlab = NULL,
                      		pfilter = NULL,
                       		shade.col = "grey80",
                       		num = 1,
                       		hide.col = NULL,
                       		lorder = NULL,
													plot = TRUE,
                       		verbose = FALSE, ...){

	##CHECK
  # M by M matrix
#	if ( dim ( d ) [ 1 ] == dim ( d ) [ 2 ] )
  # M by M^l matrix

#	if ( as.integer ( log ( dim(array, base= 2))==(log(8,base=2))


 # d <- matrix ( d, nrow = ( dim ( d )[ 1 ] ), ncol = ( dim ( d )[ 2 ] ))

	
	d <- as.matrix (d) 
	
  trmatplot.default ( d = d, seed = seed, 
											rowconstraint = rowconstraint, morder = morder,
											cspal = cspal, cpal = cpal, main = main, 
											xlab = xlab, ylab = ylab, ylim = ylim, 
											xtlab = xtlab, ytlab = ytlab,
                      pfilter = pfilter,
                      shade.col = shade.col,
                      num = num,
                      hide.col = hide.col,
                      lorder = lorder,
											plot = plot,
                      verbose = verbose, title, ... )

}
