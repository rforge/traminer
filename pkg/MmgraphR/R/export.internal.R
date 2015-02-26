##' -------------------------------------------------------- #
##' Author:          Pauline (Poulcheria) Adamopoulou, padamopo@gmail.com
##' Date:            2014-12-03
##'
##' Description:
##' Internal functions exported by main functions
##'
##'
##'	Contents:
##'	depmix.fitted.trmat	: extract the  probability transition matrix
##'												from an object of class 'depmix.fitted'
##'	ismax	:					is maximum
##'	ismin	:					is minimum
##'	paths	:					create paths to be used
##'	smax	:					state maximum
##'	smin	:					state minimum
##'
##' sum1	:					check if sum of each row of matrix equals 1
##'	tmax	:					total maximum 
##'	tmin	:					total minimum
##' 
##' -------------------------------------------------------- #

##
# depmix.fitted.trmat
##
trmat.depmix.fitted <- function ( d ) {
	
  M <- attributes ( d ) $ nstates

  Mat <- matrix ( 0, M, M )
	
  for ( i in 1 : M ) {
		
    for ( j in 1 : M ) {
		
      Mat [ i, j ] <- ( attributes ( fm ) $ transition [[ i ]] )@ parameters $ coefficients [ j ]
			
    } 

  } 
  
  Mat

}

##
# ismax
##
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


##
# ismin
##
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


##
# paths
##
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


##
# smax
##
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

##
# smin
##
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

##
# sum1
##
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

##
# tmax
##
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

##
# tmin
##
tmin <- function ( M , l , dt , shade.col , num ){
  
  ml <- M ^ l
  
  dt <- dt [ order ( dt $ w , decreasing = FALSE ) , ]
  
  vecdatch <- as.vector ( dt$ch )
  
  for ( i in ( num +1 ) : ( ml * M ) ){
    
    vecdatch [ i ] <- shade.col
    
  }
  
  ch <- vecdatch
  
  dt$ch <- as.data.frame ( ch )  
  
  dt
}
