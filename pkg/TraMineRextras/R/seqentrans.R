## Adds 2 variables to the data element of a seqelist object
##   ntrans : the length of the event sequence (number of transitions)
##   nevent : the total number of events

seqentrans <- function(fsubseq){
     fsubseq$data$ntrans <- sapply(strsplit(as.character(fsubseq$subseq), "-"),
     length)
     fsubseq$data$nevent <- fsubseq$data$ntrans - 1 +
     sapply(strsplit(as.character(fsubseq$subseq), ","), length)
     return(fsubseq)
}
