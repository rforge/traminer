## returns the modal value of a variable
## author: Gilbert Ritschard

rowmode <- function(v, except = "*") {
	tt <- table(v)
	tt <- tt[!(names(tt) %in% except)]
	return(names(tt)[which.max(tt)])
}
