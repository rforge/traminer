
debuglevel <- function(level=NULL) {
	if(is.null(level)){
		return(.Call(SD_getSDDebugLevel))
	}
	.Call(SD_setSDDebugLevel,as.integer(level))
	return(level)
	
}
