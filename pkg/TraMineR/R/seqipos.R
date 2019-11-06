## proportion of positive spells (or states)

seqipos <- function(seqdata, dss=TRUE, pos.states=NULL, with.missing=FALSE){

	if (!inherits(seqdata,"stslist"))
		stop(" [!] data is NOT a sequence object, see seqdef function to create one")
  if (is.null(pos.states))
		stop(" [!] pos.states is null, it should contain at least one valid state")

  alph <- alphabet(seqdata)
  void <- attr(seqdata,"void")
	nr <- attr(seqdata,"nr")
  if (with.missing)
    alph <- c(alph,nr)
  if (!all(pos.states %in% alph)){
    stop(" [!] invalid values in pos.states: ",paste(pos.states[!pos.states %in% alph], collapse=",
    "))
  }

  if (length(pos.states)!=length(unique(pos.states)))
    stop(" [!] Multiple occurrences of same state in pos.states")

  if (with.missing)
    recodes <- list("p"=pos.states, "!#%" = void)
  else
    recodes <- list("p"=pos.states, "*" = c(nr,void))

  if (dss)
    s <- seqdss(seqdata, with.missing = with.missing)
  else
    s <- seqdata

  sbinary <- suppressWarnings(seqrecode(s, recodes = recodes, otherwise="m"))
  npos <- rowSums(sbinary=="p")
  nneg <- rowSums(sbinary=="m")
  ratio <- npos/(nneg + npos)
  return(ratio)
}
