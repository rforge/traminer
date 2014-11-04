#'@keywords internal
smax <- function ( M , l , dt , shade.col ){
  
  ml <- M ^ l
  
  vecdatch <- as.vector ( dt$ch )
  
  for ( i in 1 : ( ml * M ) ){
    
    if ( dt $ mx [ i ] == "FALSE" ) {
      
      vecdatch [ i ] <- shade.col 
    }
    else {
      vecdatch [ i ] <- as.vector (dt $ ch )[ i ] 
    }
  }
  ch <- vecdatch
  
  dt$ch <- as.data.frame ( ch )  
  
  dt
}
