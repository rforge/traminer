#'@keywords internal
trmatplot.depmix.fitted <- function(d, seed = NULL, 
													rowconstraint = TRUE, morder = 1,
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

	## ----------------------CHECK
  ## probability transition matrix within depmix.fitted object
	if (!class(d@transition[[1]]@parameters$coefficients)=="numeric") {
		stop("[!] ensure that the depmix.fitted object contains a probability transition matrix")
	}
	## ----------------------PREPARE
	d <- trmat.depmix.fitted(d)
 
	trmatplot.default(d = d, seed = seed, 
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
                    verbose = verbose, ...)
}

