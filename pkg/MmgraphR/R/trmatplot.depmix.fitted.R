#'@keywords internal
trmatplot.depmix.fitted <- function (d, seed = NULL, cspal = NULL, cpal = NULL, title = NULL,
                	       	xlab =  NULL, ylab = NULL, ylim = NULL, xtlab = NULL,
                      		pfilter = NULL,
                       		shade.col = "grey80",
                       		num = 1,
                       		hide.col = NULL,
                       		lorder = NULL,
                       		verbose = FALSE, ...){

  d <- depmix.fitted.trmat ( d )
 
  trmatplot.default ( d = d, seed = seed, cpal = cpal, cspal = cspal, title = title, 
			xlab = xlab, ylab = ylab, ylim = ylim, xtlab = xtlab,
                        pfilter = pfilter,
                        shade.col = shade.col,
                        num = num,
                        hide.col = hide.col,
                        lorder = lorder,
                        verbose = verbose, ... ) 
  
}
