##' -------------------------------------------------------- #
##' Author:          Pauline (Poulcheria) Adamopoulou, padamopo@gmail.com
##' Date:            2014-11-24
##'
##' Description:
##' Unexported functions imported from other packages
##'
##' References:
##' march:      https://r-forge.r-project.org/projects/march/
##'
##' Contents:
##' march:                       march.dcmm.h.compactA
##' -------------------------------------------------------- #

##
##
march.dcmm.h.compactA <- function ( d ){
  RA <- array ( 0, c ( d@M ^ d@orderHC, d@M ) )
  
  for ( r2 in 1 : d@M ) {
    for ( r1 in 1 : d@M ^ (d@orderHC-1) ){
      for ( c in 1:d@M ){
        RA [ ( r1 - 1 ) * d@M + r2, c ] <- d@A [ ( r1 - 1 ) * d@M + r2, r1 + (c-1) * d@M ^ ( d@orderHC - 1 ) ]
      }
    }
  }
  
  RA
}
#author: Ogier Maitre, function in march.dcmm.R lines 62-74
