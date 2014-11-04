##########################################
# Plotting a probability transition matrix
##########################################

trmat<-matrix( c (0.1, 0.05, 0.05, 0.80,
                  0.06, 0.02, 0.03, 0.89,
                  0.03, 0.01, 0.01, 0.95,
                  0, 0, 0, 1), nrow = 4, ncol = 4, byrow = TRUE)

trmatplot(trmat)

#--- Setting a seed so that the graphic can be replicated

trmatplot(trmat, seed = 2)

###########################################################
# pfilter: Filtering out most (or least) probable sequences
###########################################################

#--- The most probable sequence given a state

trmatplot(trmat, seed = 2, pfilter = "smax")

#--- The least probable sequence given a state

trmatplot(trmat, seed = 2, pfilter = "smin")

#--- The four most probable sequnces

trmatplot(trmat, seed = 2, pfilter = "tmax", num = 4 )

#--- The ten least probable sequences

trmatplot(trmat, seed = 2, pfilter = "tmin", num = 10 )

######################################################
# filter: Highlighting a specific sequence of interest
######################################################

trmatplot(trmat, seed = 2, filter = list(type = "sequence", value = "(1)-(4)"))
