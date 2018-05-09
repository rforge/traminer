##'--------------------------------------------------------- #
##' Author:          Pauline (Poulcheria) Adamopoulou
##' E-Mail:          padamopo@gmail.com
##' Date:            2018-05-08
##'
##' Description:				New generics for exported methods
##' 
##'
##' Overview:
##' trmatplot:          parallel coordinate plot for a
##'											probability transition matrix		
##'
##' Last modifications:
##' 2018-05-08: include argument clustername for objects of class mhmm
##' 2018-05-08: include argument main because argument title is deprecated
##' 2018-05-02: include argument morder in 'trmatplot'
##' 2014-12-05: include all arguments in 'trmatplot' 
##' 2014-11-07: added 'trmatplot' generic
##'--------------------------------------------------------- #

trmatplot <- function( d, seed = NULL, rowconstraint = TRUE, morder = 1,
											cspal = NULL, cpal = NULL, main = NULL,
                      xlab =  NULL, ylab = NULL, ylim = NULL, xtlab = NULL, ytlab = NULL,
											pfilter = NULL, shade.col = "grey80", num = 1,
                      hide.col = NULL, lorder = NULL, plot = TRUE, verbose = FALSE, ...) UseMethod("trmatplot")

##
# depmix.fitted.trmat
##
trmat.depmix.fitted <- function ( d ) UseMethod("trmat")
