#'@keywords internal
trmatplot.msm <- function(d, seed = NULL, 
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

	## ----------------------PREPARE
	d <- pmatrix.msm(d)

	## ----------------------
	if (verbose) {
		cat(" [>] extracting the probability transition matrix to be plotted using the default settings of the function msm::pmatrix.msm\n")
	}
 
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

