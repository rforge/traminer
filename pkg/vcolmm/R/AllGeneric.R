## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2014-03-17
##
## Description:
## Generics for methods for S3 and S4 classes
##
## Overview:
## S3 Methods for olmm class:
## neglogLik       negative log Likelihood
## ranefCov        covariance matrix for random effects
## --------------------------------------------------------- #

## S3 generics

neglogLik <- function(object, ...) UseMethod("neglogLik")

ranefCov <- function(object, ...) UseMethod("ranefCov")

splitpath <- function(tree, ...) UseMethod("splitpath")

## S4 generics

