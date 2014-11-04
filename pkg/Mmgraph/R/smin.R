#'@keywords internal
smin <- function ( M , l , dt , shade.col ){
  
  ml <- M ^ l
  
  vecdatch <- as.vector ( dt$ch )
  
  for ( i in 1 : ( ml * M ) ){
    
    if ( dt $ mn [ i ] == "TRUE" ) {
      
      vecdatch [ i ] <- as.vector (dt $ ch )[ i ] 
    }
    else {
      vecdatch [ i ] <- shade.col
    }
  }
  ch <- vecdatch
  
  dt$ch <- as.data.frame ( ch )  
  
  dt
}
