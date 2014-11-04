#'@keywords internal
ismin <- function (M, l, dt) {
  ml <- M^l
  k <- 0 
  for ( j in ( 0 : (ml-1) ) ){
    
    for ( i in  ( 1 : M ) ){
      
      if ( dt $ w [ i + k ] == min ( dt $ w [ (j * M + 1) : (( j + 1) * M) ]) ) {
        dt $ mn [ i + k ] <- "TRUE"
      } 
      
      else {
        dt $ mn [ i + k ] <- "FALSE" 
      }
    }
    k <- k + M
  }
  dt
}
