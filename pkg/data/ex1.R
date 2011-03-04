## Example data file with missing values

s1 <- c(NA,NA,NA,"A","A","A","A","A","A","A","A","A","A")
s2 <- c("D","D","D","B","B","B","B","B","B","B",NA,NA,NA)
s3 <- c(NA,"D","D","D","D","D","D","D","D","D","D",NA,NA)
s4 <- c("A","A",NA,NA,"B","B","B","B","D","D",NA,NA,NA)
s5 <- c("A",NA,"A","A","A","A",NA,"A","A","A",NA,NA,NA)
s6 <- c(NA,NA,NA,"C","C","C","C","C","C","C",NA,NA,NA)
s7 <- rep(NA,13)

ex1 <- rbind(s1,s2,s3,s4,s5,s6,s7)

w <- c(17.7, 4.9, 5.1, 29.3, 0.7, 2.3, 0)
ex1 <- data.frame(ex1, w)

names(ex1) <- c(paste("[P",1:13,"]",sep=""), "weights")

