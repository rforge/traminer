#'@keywords internal
trmatplot.array <- function(d, seed = NULL, 
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

	## ----------------------CHECK
  ## M by M^l matrix
	if (ncol(d) ^ morder != length(as.numeric(d)) / ncol(d)) {
		stop("[!] check that 'morder' (the order of the probability transition matrix) is correctly specified")
	}

	## ----------------------PREPARE
	d <- matrix(as.numeric(d), nrow = ncol(d)^morder, ncol = ncol(d), byrow = TRUE)
	
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
