## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2014-05-02
##
## Description:
## Generics for methods for S3 and S4 classes
##
## Overview:
## S3 Methods for olmm:
## ranefCov:        extracts covariance matrixof models with
##                  random effects.
## oobrisk:         estimates out-of-bag risk.
## otsplot:         ordinal time series plot.
## splitpath:       extracts the splitting path of a tree
##                  structure.
## stabpath:        computes stability paths for variable
##                  selection.
## --------------------------------------------------------- #

## S3 generics

ranefCov <- function(object, ...) UseMethod("ranefCov")

extract <- function(object, ...) UseMethod("extract")

cvrisk <- function(object, ...) UseMethod("cvrisk")

oobrisk <- function(object, ...) UseMethod("oobrisk")

otsplot <- function(x, ...) UseMethod("otsplot")

splitpath <- function(tree, ...) UseMethod("splitpath")

stabpath <- function(object, q, ...) UseMethod("stabpath")

