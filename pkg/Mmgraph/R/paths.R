#'@keywords internal
paths <- function ( M , l ) {
  ml <- M ^ l
  
  pathsmatrix <- matrix ( 0 , ml*M , l+1)
  
  pathsmatrix [ , l + 1 ] <- rep ( 1 : M , ml )
  
  if ( l >= 1 ) {
    for ( j in 1 : l ){
      pathsmatrix [ , j ] <- rep ( 1 : M , M ^ (l-j) , each = M ^ j)
    }
    
  }
  as.data.frame ( pathsmatrix )
}
