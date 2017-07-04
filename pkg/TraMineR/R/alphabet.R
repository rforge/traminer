## ============================================
## Retrieve the alphabet from a sequence object
## ============================================

alphabet <- function(seqdata) {

	if (inherits(seqdata,"stslist")){
    statl <- attr(seqdata,"alphabet")
  }
  else if (inherits(seqdata,"seqelist")){
    statl <- levels(seqdata)
  }
  else {
		stop("seqdata is nor a state sequence object, nor en event sequence object. Use seqdef or seqecreate.")
  }

return(statl)
}
	
