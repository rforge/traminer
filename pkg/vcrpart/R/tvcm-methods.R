## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2014-05-05
##
## Description:
## S3 methods for tvcm objects
##
## References:
## party:           http://CRAN.R-project.org/package=party
## partykit:        http://CRAN.R-project.org/package=partykit
##
## Methods:
## coef, coefficients:
## extract:
## fitted:
## formula:
## getCall:             extract original call
## logLik:              extract log Likelihood
## model.frame:
## nobs:
## predict:             predict responses (see prediction of 'olmm' class)
## print:               print tvcm objects
## prune:
## ranef:
## resid, residuals:    extract residuals
## splitpath:
## weights:
## --------------------------------------------------------- #

coef.tvcm <- function(object, ...) tvcm_get_estimates(object, ...)

coefficients.tvcm <- coef.tvcm 

deviance.tvcm <- function(object, ...)
  return(deviance(extract(object, "model")))

extract.tvcm <- function(object, what = c("control", "model", 
                                   "sctest", "p.value",
                                   "riskgrid", "selected", 
                                   "coef", "sd", "var"),
                         step = NULL, ...) {
  
  what <- match.arg(what)
  splitpath <- object$info$splitpath
  if (length(splitpath) > 0 && is.null(step))
    step <- 1:length(splitpath)
  sctest <- object$info$control$method == "mob"
  rval <- NULL
  
  if (what == "control") {

    return(object$info$control)
    
  } else if (what == "sctest" && sctest) {
    
    rval <- lapply(splitpath[step], function(x) x$sctest)
    return(rval)
    
  } else if (what == "riskgrid" && !is.null(splitpath)) {
    
    rval <- lapply(splitpath[step], function(x) x$riskgrid)
    return(rval)
    
  } else if (what == "model") {
    
    return(object$info$model)
    
  } else if (what == "selected" && depth(object) > 0) {

    node <- setdiff(nodeids(object), nodeids(object, terminal = TRUE))
    varids <- unlist(nodeapply(object, node, function(node) node$split$varid))
    if (length(varids) > 0L)
        rval <- unique(colnames(object$data)[varids])
    return(rval)
    
  } else if (sctest && what == "p.value" && !is.null(splitpath)) {
    
    rval <- unlist(sapply(splitpath[step],
                          function(x) {
                            if (is.null(x$sctest)) return(NA)
                            rval <- na.omit(c(x$sctest$p.value))
                            if (length(rval) == 0) return(NA)
                            return(min(rval, na.rm = TRUE))
                          }))
    return(rval)
    
  } else if (what %in% c("coef", "sd", "var")){

    return(tvcm_get_estimates(object, what = what))   
  }
  return(rval)
}

fitted.tvcm <- function(object, ...) {
  args <- append(list(object = object), list(...))
  args$newdata <- NULL # delete the newdata argument 
  return(do.call(predict, args = args)) # ... and call predict
}

formula.tvcm <- function(x, which = c("original", "root", "tree"), ...) {
  which <- match.arg(which)
  return(formula(x$info$formula[[which]]))
}

getCall.tvcm <- function(x, ...) return(x$info$call)

logLik.tvcm <- function(object, ...)
  return(logLik(extract(object, "model")))

model.frame.tvcm <- function(formula, ...) {
    rval <- cbind(model.frame(formula$info$model), formula$data)
    attr(rval, "terms") <- attr(formula$data, "terms")
    attr(rval, "na.action") <- attr(formula$data, "na.action")
    return(rval)
}

nobs.tvcm <- function(object, ...) nobs(extract(object, "model"), ...)

predict.tvcm <- function(object, newdata = NULL,
                         type = c("link", "response", "prob", "class",
                           "node", "coef", "ranef"),
                         ranef = FALSE, na.action = na.pass, ...) {

  ## match type
  type <- match.arg(type)

  ## resolve conflicts with the 'ranef' argument
  if (!is.null(newdata) && is.logical(ranef) && ranef)
    stop("'ranef' should be 'FALSE' or a 'matrix' if 'newdata' is not 'NULL'.")
  if (type == "ranef" & (!is.logical(ranef) | is.logical(ranef) && ranef))
    stop("for 'type = 'ranef'' the argument 'ranef' must be 'FALSE'.")
  if (type == "ranef" & !is.null(newdata))
    stop("prediction for random effect for 'newdata' is not implemented.")
  
  if (type == "ranef") return(ranef(object$info$model, ...))
  
  ## the terminal node identifiers
  ids <- nodeids(object, terminal = TRUE)
  
  ## extract node ids
  if (is.null(newdata)) {
    fitted <- object$fitted[["(fitted)"]]
    names(fitted) <- rownames(object$data)
  } else {
    fitted <- newdata$Node <- tvcm_get_node(object, newdata)
  }

  ## return fitted node ids
  if (type == "node") {

    ## return fitted nodes
    return(fitted)
  } else if (type == "coef") {
      
    ## extract individual effects
    what <- list(...)$what
    if (is.null(what)) what <- "coef"
    coef <- extract(object, what)
    fe <- vc <- re <- NULL
    
    if (!is.null(coef$fe))
      fe <- matrix(coef$fe, nrow = length(fitted),
                   ncol = length(coef$fe), byrow = TRUE,
                   dimnames = list(names(fitted),
                     names(coef$fe)))

    if (!is.null(coef$vc)) {
      vc <- lapply(fitted, function(x) coef$vc[as.character(x), ])
      vc <- matrix(unlist(vc), nrow = length(fitted), byrow = TRUE,
                   dimnames = list(names(fitted),
                     colnames(coef$vc)))
    }
    
    if (!is.null(coef$re))
      re <- matrix(coef$re, nrow = length(fitted),
                   ncol = length(coef$re), byrow = TRUE,
                   dimnames = list(names(fitted),
                     names(coef$re)))
    
    
    rval <- cbind(fe, vc, re)
    rownames(rval) <- names(fitted)
    
    return(rval)
  } else {

    ## call predict.olmm
    return(predict(object$info$model, newdata = newdata, type = type,
                   na.action = na.action, ...))
  }
  return(fitted)
}

tvcm_print <- function(x, type = c("print", "summary"), ...) {

  type <- match.arg(type)
  coef <- extract(x, "coef")
  sd <- extract(x, "sd")
  
  header_panel <- function(x) {
    rval <- x$info$title
    rval <- c(rval, "", paste("  Family: ", x$info$family$family,
                              " ", x$info$family$link, sep = ""))
    if (length(x$info$formula$original) > 0L) {
      formula <- deparse(x$info$formula$original)
      formula[1L] <- paste(" Formula: ", formula[1L], sep = "")
      rval <- c(rval, formula)
    }
    if (length(str <- deparseCall(x$info$call$data)) > 0L)
      rval <- c(rval, paste("    Data: ", str, sep = ""))
    if (length(str <- deparseCall(x$info$call$subset)) > 0L)
      rval <- c(rval, paste("  Subset: ", str, sep = ""))
    meth <- paste("  Method: ", x$info$control$method, sep = "")
    if (x$info$control$method == "mob") {
      meth <- paste(meth, " (alpha = ", format(x$info$control$alpha, ...), sep = "")
      if (x$info$control$bonferroni) 
        meth <- paste(meth, ", Bonferroni corrected", sep = "")
    } else {
      meth <-
        paste(meth, " (maxwidth = ", format(x$info$control$maxwidth, ...), sep = "")
    }                    
    meth <- paste(meth, ", minbucket = ", x$info$control$minbucket, ")", sep = "")
    rval <- c(rval, meth, "")
      
    if (type == "summary") {
      lLik <- logLik(x)
      AICtab <- cbind(AIC = AIC(lLik),
                      BIC = BIC(lLik),
                      logLik = as.vector(lLik),
                      deviance = deviance(x))
      rval <- c(rval, "Goodness of fit:")
      rval <- c(rval, unlist(strsplit(formatMatrix(AICtab, ...), "\n")), "")
    } 
    
    if (length(coef$re) > 0L) {
      rval <- c(rval, "Random effects:")
      VarCorr <- VarCorr(extract(x, "model"))
      rval <- c(rval, attr(VarCorr, "title"))
      rval <- c(rval, unlist(strsplit(formatMatrix(VarCorr, ...), "\n")), "")
    }
      
    if (length(coef$fe) > 0L) {
      rval <- c(rval, "Fixed effects:")
      if (type == "print") {
        coefMat <- matrix(coef$fe, 1)
        colnames(coefMat) <- names(coef$fe)
      } else {
        coefMat <- cbind("Estimate" = coef$fe,
                         "Std. Error" = sd$fe,
                         "t value" = coef$fe / sd$fe)
      }
      rval <- c(rval, unlist(strsplit(formatMatrix(coefMat, ...), "\n")), "")
    }
      
    rval <- c(rval, "Varying effects:")
    if (depth(x) ==0L && length(coef$vc) > 0L) {
      if (type == "print") {
              coefMat <- matrix(coef$vc["1", ], 1)
              colnames(coefMat) <- colnames(coef$vc)
            } else {
              coefMat <- cbind("Estimate" = coef$vc["1",],
                               "Std. Error" = sd$vc["1",],
                               "t value" = coef$vc["1",] / sd$vc["1",])
            }
      rval <- c(rval, unlist(strsplit(formatMatrix(coefMat, ...), "\n")))
    }
    rval <- c(rval, "")
    return(rval)
  }
  
  terminal_panel <- function(node) {

    nid <- as.character(id_node(node))
    if (type == "print") {
      coefMat <- matrix(coef$vc[nid, ], 1)
      colnames(coefMat) <- colnames(coef$vc)
    } else {
      coefMat <- cbind("Estimate" = coef$vc[nid, ],
                       "Std. Error" = sd$vc[nid, ],
                       "t value" = coef$vc[nid, ] / sd$vc[nid, ])
    }
    return(c("", unlist(strsplit(formatMatrix(coefMat, ...), "\n"))))
  }

  footer_panel <- function(x) {
      rval <- ""
      if (nzchar(mess <- naprint(attr(x$data, "na.action")))) 
          rval <- c(rval, paste("(", mess, ")", sep = ""), "")
      if (length(rval) == 1 && rval == "") rval <- character()
      return(rval)
  }
  
  print.party(x, header_panel = header_panel,
              terminal_panel = terminal_panel,
              footer_panel = footer_panel, ...)
}

print.tvcm <- function(x, ...)
  tvcm_print(x, type = "print", ...)

summary.tvcm <- function(object, ...)
  tvcm_print(object, type = "summary", ...)

prune.tvcm <- function(tree, alpha = NULL,
                       maxdepth = NULL, maxwidth = NULL,
                       minsplit = NULL, minbucket = NULL, 
                       nselect = NULL, maxstep = NULL,
                       terminal = NULL, ...) {

  keepids <- list(...)$keepids
  if (is.null(keepids)) keepids <- FALSE

  args <- list()
  args$data <- model.frame(tree)
  args$family <- tree$info$family
  args$weights <- tree$fitted[,"(weights)"]
  args <- append(args, tree$info$dotargs)
  
  if (!is.null(alpha) && tree$info$control$method == "partreg") {
    warning("'alpha' is not a tuning parameter for method 'partreg'")
    alpha <- NULL
  }
  
  ## prune the tree structure
  node <- tvcm_prune_node(tree, keepids, alpha, maxdepth, maxwidth,
                             minsplit, minbucket,
                             nselect, maxstep, terminal)
  
  ## if something changes ...
  if (!identical(node, tree$node)) {
    tree$node <- node
    if (depth(tree$node) > 0L) {
        args$formula <- tree$info$formula$tree
        tree$fitted[, "(fitted)"] <-
            fitted_node(tree$node, data = tree$data)
        args$data$Node <-
            tvcm_get_node(tree, tree$data, TRUE, tree$fitted[,"(weights)"])
        args$contrasts <- appendDefArgs(list(Node = contrasts(args$data$Node)),
                                      args$contrasts)
    } else {
        args$formula <- tree$info$formula$root
    }    
    tree$info$model <-
      suppressWarnings(try(do.call(tree$info$fit, args), silent = TRUE))
    if (inherits(tree$info$model, "try-error"))
      stop("tree model fitting failed")      
    if (!is.null(terminal)) {
      tree$info$splitpath <- list()
      class(tree$info$splitpath) <- "splitpath.tvcm"
    } else {
      nstep <- width(tree) - 1L
      tree$info$nstep <- nstep
      splitpath <- tree$info$splitpath
      splitpath <- splitpath[1:(nstep + 1L)]
      splitpath[[nstep+1L]]$varid <- NULL
      splitpath[[nstep+1L]]$nodeid <- NULL
      splitpath[[nstep+1L]]$cutid <- NULL
      splitpath[[nstep+1L]]$logLik <- logLik(tree$info$model)
      splitpath[[nstep+1L]]$risk <- tree$info$control$riskfun(tree$info$model)
      splitpath[[nstep+1L]]$riskgrid <- NULL
      tree$info$splitpath <-
        tvcm_update_splitpath(splitpath, tree$node, tree$data, tree$info$control)
    }
  }
  if (!is.null(alpha)) 
    tree$info$control$alpha <- alpha
  if (!is.null(maxdepth))
    tree$info$control$maxdepth <- maxdepth
  if (!is.null(maxwidth))
    tree$info$control$maxwidth <- maxwidth
  if (!is.null(minsplit))
    tree$info$control$minsplit <- minsplit
  return(tree)
}

ranef.tvcm <- function(object, ...)
  return(ranef(object$info$model, ...))

resid.tvcm <- function(object, ...)
  return(resid(object = object$info$model, ...))

residuals.tvcm <- resid.tvcm

splitpath.tvcm <- function(tree, step = 1L:tree$info$nstep, ...) {
  rval <- tree$info$splitpath[step]
  class(rval) <- "splitpath.tvcm"
  return(rval)
}
  
print.splitpath.tvcm <- function(x, ...) {  
  for (i in 1L:length(x)) {
    if (i != 1L) cat("\n")
    cat("Step:", x[[i]]$step)
    if (is.null(unlist(x[[i]]$varid))) {
      cat(" (no splitting processed)\n")
    } else {
      cat("\nSplitting variable:", x[[i]]$var)
      cat("\nNode:", x[[i]]$nodeid)
      cat("\nCutpoint: ")
      cat(paste("{",paste(x[[i]]$cutpoint, collapse = "}, {"), "}\n", sep = ""))
    }
     if (!is.null(x[[i]]$sctest) | !is.null(x[[i]]$riskgrid))
       cat("\nDetails:\n")
    if (!is.null(x[[i]]$sctest)) {
      cat("\nCoefficient constancy tests (p-value):\n")
      print(x[[i]]$sctest$p.value, ...)
    }
    if (!is.null(x[[i]]$riskgrid)) {
      cat("\nRisk-minimizing grid search:\n") 
      for (j1 in 1:length(x[[i]]$riskgrid)) {
        for (j2 in 1:length(x[[i]]$riskgrid[[j1]])) {
          cat("Variable:", names(x[[i]]$riskgrid)[j1],
              "Node:", sub("Node", "", names(x[[i]]$riskgrid[[j1]])[j2]), "\n\n")
          print(x[[i]]$riskgrid[[j1]][[j2]], ...)
        }
      }
    }
  }
}

weights.tvcm <- function(object, ...) {
  weights(extract(object, "model"))
}
