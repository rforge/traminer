##' -------------------------------------------------------- #
##' Author:          Pauline (Poulcheria) Adamopoulou, padamopo@gmail.com
##' Date:            2018-06-05
##'
##' Description:
##' Unexported functions within other packages imported here 
##'
##' References:
##' march:      https://cran.r-project.org/package=march
##'
##' Contents:
##' march:					march.dcmm.h.compactA
##'	
##' -------------------------------------------------------- #

##
##
march.dcmm.h.compactA <- function(d) {
# Computes the reduced form of the probability transition matrix of hidden states.
  #
  # Args:
  #   d: An object of class 'march.Dcmm'.
  #
  # Returns:
  #   The reduced form of the probability transition matrix of hidden states 
  #		as an object of class 'matrix'.
  RA <- array(0, c(d@M ^ d@orderHC, d@M))
  
  for (r2 in 1:d@M) {
    for (r1 in 1:d@M ^ (d@orderHC-1)) {
      for (c in 1:d@M ){
        RA[(r1 - 1) * d@M + r2, c] <- 
							d@A[(r1 - 1) * d@M + r2, r1 + (c-1) * d@M ^ (d@orderHC - 1)]
      }
    }
  }
  return(RA)
}
#author: Ogier Maitre, function in march.dcmm.R lines 62-74
