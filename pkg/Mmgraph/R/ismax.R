#'@keywords internal
ismax <- function (M, l, dt) {
  ml <- M^l
  k <- 0 
  for ( j in ( 0 : (ml-1) ) ){
    
    for ( i in  ( 1 : M ) ){
      
      if ( dt $ w [ i + k ] == max ( dt $ w [ (j * M + 1) : (( j + 1) * M) ]) ) {
        dt $ mx [ i + k ] <- "TRUE"
      } 
      
      else {
        dt $ mx [ i + k ] <- "FALSE" 
      }
    }
    k <- k + M
  }
  dt
}
