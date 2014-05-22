## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2014-05-15
##
## Description:
## General utility functions 
##
## Overview:
## appendDefArgs:          over-write default arguments.
## addEmptyChar:           add empty spaces to a character.
## formatMatrix:           format matrices for print functions.
## vcrpart_slot
## vcrpart_value_space:    extract the values space of data data.frame.
## vcrpart_formula_eta:    extracts a list of predictor formulas
## vcrpart_formula_cond:   extracts a list of predictor formulas
## vcrpart_formula:        extracts a list of predictor formulas
## renameCoefs:            modify coefficient labels
## --------------------------------------------------------- #

appendDefArgs <- function(args, default) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Over-writes default- with user-specified arguments
  ##
  ## Arguments:
  ## args:    a list of user-specified arguments.
  ## default: a list of default arguments.
  ##
  ## Value:
  ## A list with arguments.
  ## ------------------------------------------------------- #
  
  if (is.null(args)) return(default)
  subs <- setdiff(names(default), names(args))
  if (length(subs) > 0)
    for (i in subs) args[[i]] <- default[[i]] 
  return(args)
}

addEmptyChar <- function(x, nchar, justify = "left") {

  ## ------------------------------------------------------- #
  ## Description:
  ## Add empty spaces to characters to have the same length.
  ##
  ## Arguments:
  ## x:       character vector.
  ## nchar:   integer. The size per string.
  ## justify: "left" or "right". The place where the spaces
  ##          should be added.
  ##
  ## Value:
  ## a character vector.
  ## ------------------------------------------------------- #
  
  for (i in 1:length(x))
    if (justify == "left") {
      x[i] <- paste(x[i], paste(rep(" ", nchar - nchar(x[i])),
                                collapse = ""), sep = "")
    } else if (justify == "right") {
      x[i] <- paste(paste(rep(" ", nchar - nchar(x[i])),
                          collapse = ""), x[i], sep = "")
    }
  return(x)
}

deparseCall <- function(x) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Converts a call slot into a character. The function is
  ## used for summary and print functions.
  ##
  ## Arguments:
  ## x:       slot of a call, e.g., call$formula.
  ##
  ## Value:
  ## a character.
  ## ------------------------------------------------------- #
    
    if (is.null(x)) return(character())
      rval <- paste(deparse(x), collapse = "\n")
      if (grepl("structure(list(", rval, fixed = TRUE)) rval <- character()
    return(rval)
}

formatMatrix <- function(x, ...) {
  
  ## ------------------------------------------------------- #
  ## Description:
  ## Format colnames etc. of a matrix to prettify the
  ## print functions.
  ##
  ## Arguments:
  ## x: a matrix
  ##
  ## Value:
  ## a matrix.
  ## ------------------------------------------------------- #
  
  rowNames <- rownames(x)
  colNames <- colnames(x)
  if (is.null(rowNames)) rowNames <- rep("", nrow(x))
  if (!is.null(colNames)) rowNames <- append("", rowNames)
  rowNames <- addEmptyChar(rowNames, max(nchar(rowNames)))
  rval <- format(x, ...)
  if (!is.null(colNames)) {
    for (i in 1L:ncol(rval))
      if (nchar(colNames[i]) < max(nchar(rval[i])))
        colNames[i] <- addEmptyChar(colNames[i], max(nchar(rval[i])), "right")
    for (i in 1L:ncol(rval)) 
      for (j in 1L:nrow(rval)) 
        if (nchar(rval[j, i]) < nchar(colNames[i]))
          rval[j, i] <- addEmptyChar(rval[j, i], nchar(colNames[i]), "right")
    rval <- rbind(colNames, rval)
  }
  rval <- cbind(rowNames, rval)
  rval <- apply(rval, 1, paste, collapse = "  ")
  rval <- paste(rval, "\n", sep = "")
  rval <- paste(rval, collapse = "")
  return(rval)
}

vcrpart_slot <- function(object, name) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Extract a slot from a S3 or a S4 object.
  ##
  ## Arguments:
  ## object: a fitted model of class 'glm' or 'olmm'. 
  ## name:   character. The name of the slot to access.
  ##
  ## Value:
  ## The contents of the slot 'name'.
  ## ------------------------------------------------------- #
    
  rval <- NULL
  if (isS4(object) && .hasSlot(object, name))
    rval <- slot(object, name)
  if (!isS4(object) && name %in% names(object)) 
    rval <- object[[name]]
  return(rval)
}

vcrpart_value_space <- function(data, neval = 50L) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Extract the values space of data data.frame.
  ##
  ## Arguments:
  ## data:  a data.frame.
  ## neval: the maximum number of values to be evaluated
  ##        for numeric vectors.
  ##
  ## Value:
  ## A list with values of the variables in 'data'
  ## ------------------------------------------------------- #
  
  FUN <- function(x, neval) {
    rval <- sort(unique(x))
    if (is.numeric(x) && length(rval) > neval) {
      rval <- as.double(quantile(rval, seq(0,1,length.out = neval)))
    }
    if (is.factor(x)) rval <- droplevels(rval)
    return(rval)
  }
  return(lapply(as.list(data), FUN, neval = neval))
}

fe <- function(formula, intercept = TRUE) {
  if (is.logical(intercept))
    intercept <- ifelse(intercept, "ce", "none")
  if (!intercept %in% c("none", "ce"))
    stop("'intercept' in 'fe' must be one of 'none' or 'ce'")
  mc <- match.call()
  if (missing(formula)) {
    formula <- if (intercept == "none") "-1" else "1"
  } else {
    formula <- deparse(mc$formula, 500L)
  }
  eta <- as.formula(paste("~", formula))
  return(list(eta = eta, cond = ~1, intercept = intercept, type = "fe"))
}

re <- function(formula, intercept = TRUE) {
  if (is.logical(intercept))
    intercept <- ifelse(intercept, "ge", "none")
  if (!intercept %in% c("none", "ce", "ge"))
    stop("'intercept' in 'fe' must be one of 'none', 'ce' or 'ge'")
  mc <- match.call()
  formula <- deparse(mc$formula)
  formula <- strsplit(formula, "|", fixed = TRUE)[[1L]]
  if (length(formula) < 2L) stop("no grouping factor 'subject' in 're'.")
  if (length(formula) > 2L) stop("'formula' in 're' is missspecified.")
  cond <- formula(paste("~", formula[2]))
  if (length(all.vars(cond)) > 1L)
    stop("maximum one grouping factor 'subject' is allowed in 're'")
  if (formula[1L] == "") formula[1L] <- "1"
  eta <- as.formula(paste("~", formula[1L]))
  return(list(eta = eta, cond = cond, intercept = intercept, type = "re"))
}

vc <- function(..., by, intercept = TRUE) {
  if (is.logical(intercept))
    intercept <- ifelse(intercept, "ce", "none")
  if (!intercept %in% c("none", "ce", "ge"))
    stop("'intercept' in 'fe' must be one of 'none', 'ce' or 'ge'")
  mc <- match.call()
  eta <- if (missing(by)) formula(~1) else formula(paste("~", deparse(mc$by, 500L)))
  subs <- if(is.null(names(mc))) 2:length(mc) else which(names(mc) == "")[-1L]
  if (length(subs) < 1L) stop("no effect modifiers in 'vc'.")
  cond <- "~"
  for (i in subs) {
    if (i != subs[1L]) cond <- paste(cond, "+")
    cond <- paste(cond, deparse(mc[[i]]))
  }
  cond <- formula(cond)
  return(list(eta = eta, cond = cond, intercept = intercept, type = "vc"))
}

ce <- function(formula) {
  mc <- match.call()
  formula <- deparse(mc$formula, 500L)
  formula <- formula(paste("~", formula))
  return(attr(terms(formula, keep.order = TRUE), "term.labels"))
}

ge <- function(formula) {
  mc <- match.call()
  formula <- deparse(mc$formula, 500L)
  formula <- formula(paste("~", formula))
  return(attr(terms(formula, keep.order = TRUE), "term.labels"))
}


vcrpart_formula_eta <- function(x, family, env) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Sets up a list with predictor formulas.
  ##
  ## Arguments:
  ## x:      a list of formula produced by 'fe', 're'  or 'vc'
  ## family: the model family, e.g. cumulative etc.
  ## env:    environment for evaluating the formula.
  ##
  ## Value:
  ## A list with sublists for formulas. For internal use only.
  ## ------------------------------------------------------- #

  terms <- terms(x$eta, specials = c("fe", "vc", "re", "ce", "ge"), keep.order = TRUE)
  if (length(unlist(attr(terms, "specials")[c("fe", "vc", "re")])) > 0)
    stop("'", x$type, "' term contains 'fe', 'vc' or 're' terms.")
  rval <- attr(terms, "term.labels")
  termFact <- rownames(attr(terms, "factors"))
  
  if (x$type %in% c("fe", "re", "vc") & inherits(family, "family.olmm")) {

    ## extract the terms for 'olmm' objects which may include
    ## terms width 'ge()' and 'ce'

    ## extract terms
    if (length(rval) > 0L) {

      subsCe <- rval %in% termFact[attr(terms, "specials")$ce]
      ceTerms <- unlist(lapply(rval[subsCe], function(x) eval(parse(text = x))))

      subsGe <- rval %in% termFact[attr(terms, "specials")$ge]
      geTerms <- unlist(lapply(rval[subsGe], function(x) eval(parse(text = x)))) 

      if (family$family %in% c("cumulative", "adjacent")) {
        geTerms <- c(geTerms, rval[!subsGe & !subsCe])
      } else {
        ceTerms <- c(ceTerms, rval[!subsGe & !subsCe])
      }
      rval <- list(paste(ceTerms, collapse = "+"),
                   paste(geTerms, collapse = "+"))

    } else {

      rval <- list("", "")
      
    }
    
    ## set intercepts
    if (x$intercept == "none") {
      int <- list(switch(x$type, fe = "-1", re = "-1", vc = "1"),
                  switch(x$type, fe = "1", re = "-1", vc = "1"))
    } else if (x$intercept == "ge") {
      int <- list(switch(x$type, fe = "1", re = "-1", vc = "1"),
                  switch(x$type, fe = "1", re = "1", vc = "Node"))
    } else if (x$intercept == "ce") {
      int <- list(switch(x$type, fe = "1", re = "1", vc = "Node - 1"),
                  switch(x$type, fe = "1", re = "-1", vc = "1"))
    }

    rval <- lapply(1L:2L, function(i) {
      if (rval[[i]] == "") return(int[[i]])
      if (int[[i]] != "1") return(paste(int[[i]], rval[[i]], sep = "+"))
      return(rval[[i]])
    })
    
    if (x$type == "re" && family$family == "cumulative" && rval[[1L]] != "-1")
      stop("category-specific random effects are not ",
           "available for the cumulative model.")

    if (rval[[1L]] == "-1" & rval[[2L]] == "-1")
      stop("the term '", x$type, "' is misspecified any may be dropped from ",
           "the specification")
    
    ## set formulas
    rval <- lapply(paste("~", rval, sep = "") , as.formula)

    if (x$type == "re") {
      int <- sapply(rval, function(x) attr(terms(x), "intercept"))
      if (all(int == 1L)) rval[[1L]] <- update(rval[[1L]], ~ . - 1)
    }
    
    ## incorporate 'Node'
    if (x$type == "vc") {

      ## add 'Node' to all terms
      addNode <- function(x) {
        terms <- attr(terms(x, keep.order = TRUE), "term.labels")
        termsNew <- paste("Node", terms, sep = ":")
        if (length(terms) > 0L) {
          x <- update(x, paste("~.-",paste(terms, collapse = "-")))
          x <- update(x, paste("~.+",paste(termsNew, collapse = "+")))
        }
        return(x)
      }
      rval <- lapply(rval, addNode)
    }
    
  } else {

    if (x$type == "vc") {
      terms <- attr(terms(x$eta, keep.order = TRUE), "term.labels")
      if (length(terms) > 0L) terms <- paste("Node", terms, sep = ":")
      if (x$intercept != "none") terms <- c("Node", terms)
      if (length(terms) == 0L) terms <- "1"
      rval <- list(~1, as.formula(paste("~", paste(terms, collapse = "+"))))
    } else {
      rval <- list(~1, x$eta)
    }
  }
  environment(rval[[1L]]) <- env
  environment(rval[[2L]]) <- env
  names(rval) <- c("ce", "ge")
  return(rval)
}

vcrpart_formula_cond <- function(x, family, env) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Sets up a list with predictor formulas and conditioning
  ## variables for a term 'vc', 'fe' or 're' in a formula
  ##
  ## Arguments:
  ## x:      a list of formula produced by 'fe', 're'  or 'vc'
  ## family: the model family, e.g. cumulative etc.
  ## env:    environment for evaluating the formula.
  ##
  ## Value:
  ## A list with sublists for formulas. For internal use only.
  ## ------------------------------------------------------- #
  
  ## extract the inner formula of the term
  rval <- x$cond
  environment(rval) <- env

  return(rval)
}

vcrpart_formula <- function(formula, family = cumulative(), env = parent.frame()) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Evaluates the input formula and returns a list with
  ## the necessary parts for constructing model matrices.
  ##
  ## Arguments:
  ## formula: the original formula.
  ## family:  the model family, e.g. cumulative() etc.
  ## env:     environment for evaluating the formula.
  ##
  ## Value:
  ## A list with sublists of formulas. For internal use only.
  ## ------------------------------------------------------- #

  types <- c("fe", "vc", "re")
  
  terms <- terms(formula, specials = types, keep.order = TRUE)
  if (any(sapply(attr(terms, "specials"), length) > 1L))
    stop("one of 'fe', 're' or 'vc' is duplicated in the formula")
  yName <- all.vars(formula)[1L]
  termLabs <- type <- attr(terms, "term.labels")
  termFact <- rownames(attr(terms, "factors"))
  subs <- termLabs %in% termFact[unlist(attr(terms, "specials"))]
  type[subs] <- substr(termLabs[subs], 1L, 2L)
  unknown <- !type %in% types
  
  ## check if all 'terms' use one of 'fe', 'vc' or 're'
  if (any(unknown)) {
    feTerm <- termLabs[unknown]
    if ("fe" %in% type) {
      feAdd <- eval(parse(text = attr(terms, "term.labels")[type == "fe"]))
      feTerm <- c(attr(terms(feAdd$eta, keep.order = TRUE), "term.labels"), feTerm)
      feTerm <- paste("fe(", paste(feTerm, collapse = "+"),
                      ", intercept='", feAdd$intercept,  "')", sep = "")
      termLabs[type == "fe"] <- feTerm
      termLabs <- termLabs[!unknown]
      type <- type[!unknown]
    } else {
      feTerm <- paste("fe(", paste(feTerm, collapse = "+"), ")", sep = "")
      termLabs <- termLabs[!unknown]
      type <- type[!unknown]
      termLabs <- c(termLabs, feTerm)
      type <- c(type, "fe")
    }
  }
  
  rval <- vector("list", 3L)
  names(rval) <- types
  ## extract formulas for different types of coefficients
  for (i in 1L:length(types)) {
    if (types[i] %in% type) {
      subs <- which(type == types[i])
      rval[[i]] <- eval(parse(text = termLabs[subs]))
      rval[[i]]$eta <- vcrpart_formula_eta(rval[[i]], family, env)
      rval[[i]]$cond <- vcrpart_formula_cond(rval[[i]], family, env)
    }
  }

  ## add original formula
  rval$original <- formula(formula, env = env)

  ## formula with all predictors
  getTerms <- function(x) {
    if (inherits(x, "formula")) {
      rval <- attr(terms(x, keep.order = TRUE), "term.labels")
      if ("vc" %in% type) {
        rval <- gsub("Node:", "", rval)
        rval <- rval[rval != "Node"]
      }
    } else {
      rval <- NULL
    }
    return(rval)
  }
  allTerms <- unlist(lapply(unlist(rval[1L:3L]), getTerms))
  if (length(allTerms) == 0L) allTerms <- "1"
  fAll <- paste(yName, paste(allTerms, collapse = "+"), sep = "~")
  fAll <- formula(fAll, env = env)
  rval$all <- fAll
  return(rval)
}
