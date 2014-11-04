#'@keywords internal
tmax <- function ( M , l , dt , shade.col , num ){
  
  ml <- M ^ l
  
  dt <- dt [ order ( dt $ w , decreasing = TRUE ) , ]
  
  vecdatch <- as.vector ( dt$ch )
  
  for ( i in ( num +1 ) : ( ml * M ) ){
    
    vecdatch [ i ] <- shade.col
    
  }
  
  ch <- vecdatch
  
  dt$ch <- as.data.frame ( ch )  
  
  dt
}
