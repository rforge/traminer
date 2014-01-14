## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2013-09-17
##
## Description:
## Generics for methods for S3 and S4 classes
##
## Overview:
## S4 Methods for olmm class:
## neglogLik       negative log Likelihood
## ranefCov        covariance matrix for random effects
## --------------------------------------------------------- #

## S3 generics

neglogLik <- function(object, ...) UseMethod("neglogLik")

## S4 generics

setGeneric(name = "ranefCov",
           def = function(object, ...) {standardGeneric("ranefCov")})

