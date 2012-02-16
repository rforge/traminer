## adds percent size to each factor levels
## author: Gilbert Ritschard

group.p <- function(group){
    if(!is.factor(group)){
        group <- factor(group)
    }
    ctb <- table(group)
    cprop <- prop.table(ctb)
    cprop.prct <- as.character(round(100*cprop,1))
    levels(group) <- paste(levels(group)," (",cprop.prct,"%)", sep="")
    return(group)
}

