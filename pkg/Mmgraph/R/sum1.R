#'@keywords internal
sum1 <- function ( d ) {
  
  nr <- nrow ( d )
  
  sm1 <- matrix (rep (0, nr))
  
  for ( i in 1 : nr ) {
    
    sm1 [ i ] <- sum ( d [ i , ] ) 
    
  }
  
  if ( any ( sm1 != 1 ) ){
    
    return ( FALSE )
    
  }
  
  else {
    
    return ( TRUE )
    
  }
}
