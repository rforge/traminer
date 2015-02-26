#'@keywords internal
trmatplot.depmix.fitted <- function (d, rowconstraint = TRUE, seed = NULL, 
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

  d <- trmat.depmix.fitted ( d )
 
   trmatplot.default ( d = d, rowconstraint = rowconstraint, seed = seed,
											cspal = cspal, cpal = cpal, title = title, 
											xlab = xlab, ylab = ylab, ylim = ylim, 
											xtlab = xtlab, ytlab = ytlab,
                      pfilter = pfilter,
                      shade.col = shade.col,
                      num = num,
                      hide.col = hide.col,
                      lorder = lorder,
											plot = plot,
                      verbose = verbose, ... )
  
}

