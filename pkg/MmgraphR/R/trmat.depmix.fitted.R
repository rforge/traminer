trmat <- function ( d ) {
	
  M <- attributes ( d ) $ nstates

  Mat <- matrix ( 0, M, M )
	
  for ( i in 1 : M ) {
		
    for ( j in 1 : M ) {
		
      Mat [ i, j ] <- ( attributes ( d ) $ transition [[ i ]] )@ parameters $ coefficients [ j ]
			
    } 

  } 
  
  Mat
}
