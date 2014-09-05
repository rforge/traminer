##' -------------------------------------------------------- #
##' Author:          Reto Buergin
##' E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
##' Date:            2014-09-02
##'
##' Description:
##' S3 methods for tvcm objects
##'
##' References:
##' party:           http://CRAN.R-project.org/package=party
##' partykit:        http://CRAN.R-project.org/package=partykit
##'
##' Methods:
##' coef, coefficients:
##' extract:
##' fitted:
##' formula:
##' getCall:             extract original call
##' logLik:              extract log Likelihood
##' model.frame:         extract the total model frame including model
##'                      and partitioning variables
##' nobs:                extract the number of observations
##' predict:             predict responses (see prediction of 'olmm' class)
##' print:               print tvcm objects
##' prune:               prune the tree
##' ranef:               extract the predicted random effects
##' resid, residuals:    extract residuals
##' splitpath:           show information of the splitting procedure
##' weights:             extract the weights
##'
##' Modifications:
##' 2014-09-02: added 'prunepath.tvcm' and 'print.prunepath.tvcm' functions
##' 2014-09-02: various modifications in 'prune.tvcm' for accelerations:
##'             - 'do.call' was replaced by 'eval'
##'             - new option 'keeploss' (reuses information of the previous
##'               pruning step)
##'             - new option 'papply' (even though the accrelation is not
##'               as efficient as expected)
##' 2014-08-13: modified option 'type = "coef"' of predict.tvcm
##'             to deal with multiple components
##' -------------------------------------------------------- #

coef.tvcm <- function(object, ...) tvcm_get_estimates(object, ...)

coefficients.tvcm <- coef.tvcm 

extract.tvcm <- function(object, what = c("control", "model", 
                                   "nodes", "sctest", "p.value",
                                   "lossgrid", "selected", 
                                   "coef", "sd", "var"),
                         steps = NULL, ...) {
  
  what <- match.arg(what)
  splitpath <- object$info$splitpath
  if (length(splitpath) > 0 && is.null(steps))
    steps <- seq(1L, object$info$nstep)
  rval <- NULL
  sctest <- object$info$control$sctest
  
  if (what == "control") {

    return(object$info$control)
    
  } else if (what == "sctest" && sctest) {
    
    rval <- lapply(splitpath[steps], function(x) x$sctest)
    return(rval)
    
  } else if (what == "lossgrid" && !is.null(splitpath)) {
    
    rval <- lapply(splitpath[steps], function(x) x$lossgrid)
    return(rval)
    
  } else if (what == "model") {
    
    return(object$info$model)
    
  } else if (what == "selected" && depth(object) > 0) {

    rval <-
      lapply(object$info$node,
             function(node) {
               ids <- setdiff(nodeids(node), nodeids(node, terminal = TRUE))
               rval <- unlist(nodeapply(node, ids, function(node) node$split$varid))
               if (length(rval) > 0L) rval <- unique(colnames(object$data)[rval])
               return(rval)
             })
    return(rval)
    
  } else if (sctest && what == "p.value" && !is.null(splitpath)) {
    
    rval <- unlist(sapply(splitpath[steps],
                          function(x) {
                            if (is.null(x$sctest)) return(NA)
                            rval <- na.omit(unlist(x$sctest))
                            if (length(rval) == 0) return(NA)
                            return(min(rval, na.rm = TRUE))
                          }))
    return(rval)
    
  } else if (what %in% c("coef", "sd", "var")) {

    return(tvcm_get_estimates(object, what = what, ...))   

  } else if (what == "nodes") {

    return(object$info$node)
    
  }
  return(rval)
}

extractAIC.tvcm <- function(fit, scale, k = 2, ...) {
  extractAIC(extract(fit, "model"), scale, k, ...)
}

fitted.tvcm <- function(object, ...) {
  args <- append(list(object = object), list(...))
  args$newdata <- NULL # delete the newdata argument 
  return(do.call(predict, args = args)) # ... and call predict
}

formula.tvcm <- function(x, ...) return(x$info$formula$original)


getCall.tvcm <- function(x, ...) return(x$info$call)

logLik.tvcm <- function(object, ...)
  return(logLik(extract(object, "model")))

model.frame.tvcm <- function(formula, ...) {
    rval <- cbind(model.frame(formula$info$model), formula$data)
    attr(rval, "terms") <- attr(formula$data, "terms")
    attr(rval, "na.action") <- attr(formula$data, "na.action")
    return(rval)
}

neglogLik2.tvcm <- function(object, ...)
  return(-2 * as.numeric(logLik(extract(object, "model"))))

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
  ids <- lapply(object$info$node, nodeids, terminal = TRUE)
  
  ## substitute 'newdata' by learning sample
  if (is.null(newdata)) newdata <- model.frame(object)

  ## extract node ids
  fitted <- as.data.frame(tvcm_get_node(object, newdata, formList = object$info$formula))
  newdata[, names(fitted)] <- fitted

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
    
    if (length(coef$fe) > 0L)
      fe <- matrix(coef$fe, nrow = nrow(fitted),
                   ncol = length(coef$fe), byrow = TRUE,
                   dimnames = list(rownames(fitted), names(coef$fe)))

    if (!is.null(coef$vc)) {
      for (pid in 1:length(coef$vc)) {
        if (ncol(coef$vc[[pid]]) > 0L) {
          vcPid <- lapply(fitted[, pid], function(x) coef$vc[[pid]][as.character(x), ])
          vcPid <- matrix(unlist(vcPid), nrow = nrow(fitted), byrow = TRUE)
          rownames(vcPid) <- rownames(fitted)
          colnames(vcPid) <- colnames(coef$vc[[pid]])
          vc <- cbind(vc, vcPid)
        }
      }
    }
    
    if (length(coef$re) > 0L)
      re <- matrix(coef$re, nrow = nrow(fitted),
                   ncol = length(coef$re), byrow = TRUE,
                   dimnames = list(rownames(fitted), names(coef$re)))
    
    rval <- cbind(fe, vc)
    terms <- unique(colnames(rval))
    rval <-
      sapply(terms, function(x) rowSums(rval[, colnames(rval) == x, drop = FALSE]))
    rval <- cbind(rval, re)
    rownames(rval) <- rownames(fitted)
    
    return(rval)
  } else {

    ## call predict of fitting function
    return(predict(object$info$model, newdata = newdata, type = type,
                   ranef = ranef, na.action = na.action, ...))
  }
  return(fitted)
}

tvcm_print <- function(x, type = c("print", "summary"),
                       etalab = c("int", "char", "eta"), ...) {

  type <- match.arg(type)
  etalab <- match.arg(etalab)
  coef <- extract(x, "coef")
  sd <- extract(x, "sd")
  yLevs <- if (x$info$fit == "olmm") levels(x$info$model$y) else NULL
  
  cat(x$info$title, "\n")
  cat("  Family:", x$info$family$family, x$info$family$link, "\n")
  cat(" Formula:", deparse(x$info$formula$original, 500L)[1L], "\n")  
  if (length(str <- deparseCall(x$info$call$data)) > 0L)
    cat("    Data:", str, "\n")
  if (length(str <- deparseCall(x$info$call$subset)) > 0L)
    cat("  Subset:", str, "\n")
  if (x$info$control$sctest)
    cat(paste("   Tests: alpha = ", format(x$info$control$alpha, ...),
              if (x$info$control$bonferroni) ", nodewise Bonferroni corrected\n",
              sep = ""))
  
  if (type == "summary") {
    cat("\nGoodness of fit:\n")
    lLik <- logLik(x)
    AICtab <- data.frame(AIC = AIC(lLik),
                         BIC = BIC(lLik),
                         logLik = as.vector(lLik),
                         row.names = "")
    print(AICtab, ...)
  } 
    
  if (length(coef$re) > 0L) {
    cat("\nRandom effects:\n")
    VarCorr <- VarCorr(extract(x, "model"))
    if (x$info$fit == "olmm") {
      VarCorr <- olmm_rename(VarCorr, yLevs, x$info$family, etalab)
      print.VarCorr.olmm(VarCorr, ...)
    } else {
      print(VarCorr)
    }
  }
    
  if (length(coef$fe) > 0L) {
    cat("\nFixed effects:\n")
    if (type == "print") {
      coefMat <- matrix(coef$fe, 1)
      colnames(coefMat) <- names(coef$fe) 
       
    } else {
      coefMat <- cbind("Estimate" = coef$fe,
                       "Std. Error" = sd$fe,
                       "z value" = coef$fe / sd$fe)
    }
    if (x$info$fit == "olmm")
      coefMat <- olmm_rename(coefMat, yLevs, x$info$family, etalab)
    print(coefMat, ...)
  }

  terminal_panel <- function(node, partid) {
    rval <- function(node) {
      nid <- as.character(id_node(node))
      if (type == "print") {
        coefMat <- matrix(coef$vc[[partid]][nid, ], 1)
        colnames(coefMat) <- colnames(coef$vc[[partid]])
      } else {
        coefMat <- cbind("Estimate" = coef$vc[[partid]][nid, ],
                         "Std. Error" = sd$vc[[partid]][nid, ],
                         "z value" = coef$vc[[partid]][nid, ] / sd$vc[[partid]][nid, ])
        rownames(coefMat) <- colnames(coef$vc[[partid]])
      }
      if (x$info$fit == "olmm")
        coefMat <- olmm_rename(coefMat, yLevs, x$info$family, etalab)
      return(c("", unlist(strsplit(formatMatrix(coefMat, ...), "\n"))))
    }
    return(rval)
  }
  class(terminal_panel) <- "grapcon_generator"

  vcLabs <- tvcm_print_vclabs(x)

  for (pid in seq_along(coef$vc)) {
    cat(paste("\nVarying coefficient: ", vcLabs[pid], "\n", sep = ""))
    x$node <- x$info$node[[pid]]
    print.party(x, terminal_panel = terminal_panel,
                tp_args = list(partid = pid))
    if (depth(x$info$node[[pid]]) == 0L && length(coef$vc[[pid]]) > 0L) {
      if (type == "print") {
        coefMat <- coef$vc[[pid]]
      } else {
        coefMat <- cbind("Estimate" = coef$vc[[pid]],
                         "Std. Error" = as.double(sd$vc[[pid]]),
                         "z value" = as.double(coef$vc[[pid]] / sd$vc[[pid]]))
      }
      print(coefMat, ...)
    }
  }
  
  ## print footer
  if (nzchar(mess <- naprint(attr(x$data, "na.action")))) 
    cat(paste("\n(", mess, ")\n", sep = ""))
}

print.tvcm <- function(x, ...)
  tvcm_print(x, type = "print", ...)

prune.tvcm <- function(tree, dfsplit = NULL, dfpar = 2.0,
                       direction = c("backward", "forward"),
                       alpha = NULL, maxstep = NULL, terminal = NULL,
                       papply = mclapply, keeploss = FALSE, ...) {

  ## checking arguments
  direction <- match.arg(direction)
  stopifnot(is.numeric(dfpar) && length(dfpar) == 1L)
  papplyArgs <- list(...)[names(list(...)) %in% names(formals(papply))]
  stopifnot((is.logical(keeploss) | is.numeric(keeploss)) &&
            length(keeploss) == 1L)
  keeploss <- as.numeric(keeploss)
  
  tunepar <- c(!is.null(dfsplit),
               !is.null(alpha),
               !is.null(maxstep),
               !is.null(terminal))

  if (sum(tunepar) < 1L) stop("no tuning parameter specified.")
  if (sum(tunepar) > 1L) stop("only one tuning parameter is allowed.")
  
  tunepar <- c("dfsplit", "alpha", "maxstep", "terminal")[tunepar]
  
  if (tunepar == "dfsplit") {
    stopifnot(is.numeric(dfsplit) && length(dfsplit) == 1L)

  } else if (tunepar == "alpha") {
    if (!tree$info$control$sctest) {
      warning("'alpha' is not a tuning parameter for 'tree'.")
      alpha <- NULL
    }
    stopifnot(is.numeric(alpha) && length(alpha) == 1L &&
              alpha >= 0.0 && alpha <= 1.0)

  } else if (tunepar == "maxstep") {   
    stopifnot(is.numeric(maxstep) && length(maxstep) == 1L &&
              maxstep >= 0L)
    
  } else if (tunepar == "terminal") {
    errMess <- paste("'terminal' must be a list with ",
                     length(tree$info$node), "elements.")
    if (!is.list(terminal)) {
      if (is.numeric(terminal) && length(tree$info$node) == 1L) {
        terminal <- list(terminal)
      } else {
        stop(errMess)
      }
    }
    if (is.list(terminal)) {
      if (length(terminal) != length(tree$info$node))
        stop(errMess)
      if (!all(sapply(terminal, function(x) is.null(x) | is.numeric(x))))
        stop("the elements of 'terminal' must be 'NULL' or 'numeric'")
      terminal <- lapply(terminal, as.integer)
    }
  }
  
  if (tree$info$pruned &&
      (any(tunepar %in% c("alpha", "maxstep")) |
       (tunepar == "dfsplit" & direction == "forward")))
    stop("pruning on 'alpha' or 'maxstep' is not available because ",
         "'tree' has previously been prune backwards.")

  if (tunepar %in% c("alpha", "maxstep", "terminal")) {
    
    ## prune the tree structure
    node <- tvcm_prune_node(tree, alpha, maxstep, terminal)
    
    ## refit the model if something has changed
    if (!identical(node, tree$info$node)) {
      
      ## attach nodes
      tree$info$node <- node
      tree$node <- node[[1L]]
      
      ## extract new formula
      env <- environment()
      formList <- tree$info$formula
      vcRoot <- sapply(node, width) == 1L
      ff <- tvcm_formula(formList, vcRoot, tree$info$family, env)
      
      ## overwrite node predictors
      data <- model.frame(tree)
      data[, paste("Node", LETTERS[seq_along(node)], sep = "")] <-
        tvcm_get_node(tree, tree$data, TRUE, tree$fitted[,"(weights)"], formList)
      
      ## refit model
      call <- call(name = tree$info$fit,
                   formula = quote(ff$full),
                   data = quote(data),
                   family = quote(tree$info$family),
                   weights = tree$fitted[,"(weights)"])
      for (arg in names(tree$info$dotargs)) call[[arg]] <- tree$info$dotargs[[arg]]
      tree$info$model <- suppressWarnings(try(eval(call), TRUE))
      
      if (inherits(tree$info$model, "try-error"))
        stop("tree model fitting failed.")
      
      if (tunepar == "terminal")
        tree$info$pruned <- TRUE
      
      ## update control
      tree$info$control <-
        tvcm_grow_setcontrol(tree$info$control, tree$info$model, formList, vcRoot)
    }
    
    ## update 'control'
    if (!is.null(alpha)) 
      tree$info$control$alpha <- alpha
    if (!is.null(maxstep))
      tree$info$control$maxstep <- maxstep
    
  } else if (tunepar == "dfsplit") {
    
    if (direction == "backward") {
      
      ## prune as long as there remain 'weak' links   
      run <- 1L; step <- 0L;
      cols <- c("loss", "npar", "nspl")
      prunepath <- list()
      keeplosscount <- 0
      
      while (run > 0) {
        
        run <- 0 ; step <- step + 1L;
        ncollapse <- NULL
        
        ids <- lapply(tree$info$node, function(node) {
          setdiff(nodeids(node), nodeids(node, terminal = TRUE))
        })
        
        ids00 <- lapply(seq_along(tree$info$node),
                        function(pid) {
                          unlist(nodeapply(tree$info$node[[pid]], ids[[pid]],
                                           function(node) node$info$id$original))
                        })
        
        ntab <- data.frame(
                  part = rep(seq_along(tree$info$node), sapply(ids, length)),
                  node = unlist(ids),
                  loss = rep(Inf, length(unlist(ids))),
                  npar = rep(NA, length(unlist(ids))),
                  nspl = rep(NA, length(unlist(ids))),
                  dfsplit = rep(NA, length(unlist(ids))))
              
        if (nrow(ntab) > 0L) {
          
          if (step > 1L) {
            
            ## search for the parents of the previously selecte node
            spart <- otab$part[ocollapse]
            snode <- otab$node[ocollapse]
            root <- 1L
            npath <- c(root); opath <- c(root);
            while (root != snode) {
              nkids <- unlist(nodeapply(tree$info$node[[spart]], root,
                                        function(node)
                                        sapply(node$kids, function(kids) kids$id)))
              okids <- unlist(nodeapply(tree$info$node[[spart]], nkids,
                                        function(node) node$info$id$last))
              kidkids <- nodeapply(tree$info$node[[spart]], nkids, nodeids)
              kid <- sapply(kidkids, function(x) snode %in% x)
              root <- nkids[kid]
              if (root != snode) {
                npath <- c(npath, nkids[kid])
                opath <- c(opath, okids[kid])
              }
            }
            
            ## insert results of parent models
            ntab[ntab$node %in% npath & ntab$part == spart, cols] <-
              otab[otab$node %in% opath & otab$part == spart, cols]
            

            if (keeplosscount < keeploss) {
            
              ## approximate the remaining ids from the selected partition
              otab$loss <- otab$loss + (otab$loss[ocollapse] - loss0)
              otab$npar <- otab$npar - (npar0 - otab$npar[ocollapse])
              otab$nspl <- otab$nspl - (nspl0 - otab$nspl[ocollapse])
              nids <- nodeids(tree$info$node[[spart]])
              oids <- unlist(nodeapply(tree$info$node[[spart]],
                                       nodeids(tree$info$node[[spart]]),
                                       function(node) node$info$id$last))
              subs <- ntab$node %in% nids &
                ntab$part %in% spart & is.infinite(ntab$loss)
              subs <- nids %in% ntab$node[subs]
              nids <- nids[subs]
              oids <- oids[subs]
              ntab[ntab$node %in% nids & ntab$part == spart, cols] <-
                otab[otab$node %in% oids & otab$part == spart, cols]
              keeplosscount <- keeplosscount + 1
              
            } else keeplosscount <- 0
          }
          
          subs <- which(is.infinite(ntab$loss))
          if (length(subs) > 0L) {
            
            ## reestimate all models which collapse each on inner node
            prStat <- function(i) {
              term <- lapply(seq_along(tree$info$node),
                             function(pid) {
                               if (pid == ntab$part[i]) return(ntab$node[i])
                               return(NULL)
                             })
              prTree <- try(prune(tree, terminal = term), TRUE)
              if (!inherits(prTree, "try-error"))
                return(c(prTree$info$control$lossfun(prTree),
                         extractAIC(prTree$info$model)[1L],
                         sum(sapply(prTree$info$node, width) - 1L)))
              return(c(Inf, NA, NA))
            }
            stat <- do.call(papply, append(list(X = subs, FUN = prStat), papplyArgs))
            ntab[subs, cols] <- t(sapply(stat, function(x) x))
          }
          
          if (any(ntab$loss < Inf)) {
            
            ## minimum dfsplit such that the smaller model improves the fit
            loss0 <- tree$info$control$lossfun(tree)
            npar0 <- extractAIC(tree)[1L]
            nspl0 <- sum(sapply(tree$info$node, width) - 1L)
            ntab$dfsplit <- (ntab$loss + dfpar * ntab$npar - loss0 - dfpar * npar0) /
              (nspl0 - ntab$nspl)
            
            ## prune selected inner node
            if (any(ntab$dfsplit <= dfsplit)) {
              ncollapse <- which(!is.na(ntab$dfsplit) &
                                 !is.nan(ntab$dfsplit) &
                                 ntab$dfsplit == min(ntab$dfsplit))
              if (length(ncollapse) > 1L) ncollapse <- sample(ncollapse, 1L)
              if (length(ncollapse) > 0L) {
                term <- lapply(seq_along(tree$info$node),
                               function(pid) {
                                 if (pid == ntab$part[ncollapse])
                                   return(ntab$node[ncollapse])
                                 return(NULL)
                               })
                tree <- prune(tree, terminal = term, keeploss = TRUE)
                run <- 1
              }
            }
            otab <- ntab
            ocollapse <- ncollapse
            
          } else {
            stop("fitting of nested models failed.")
          }
        }
        
        ## create a 'data.frame' for the 'prunepath' method
        tab <- data.frame(
                 part = NA,
                 node = NA,
                 loss = tree$info$control$lossfun(tree),
                 npar = extractAIC(tree)[1L],
                 nspl = sum(sapply(tree$info$node, width) - 1L),
                 dfsplit = NA,
                 row.names = "<none>")
        if (nrow(ntab) > 0L) {
          ntab$node <- unlist(ids00)
          rownames(ntab) <- seq_along(rownames(ntab))
          tab <- rbind(tab, ntab)
        }
        prunepath[[step]] <- list(step = step, tab = tab)
      }
    }
    tree$info$prunepath <- prunepath
    
  } else if (direction == "forward") {
    
    ## table with information from the latest step
    tab <- data.frame(
             step = seq(0, tree$info$nstep, 1L),
             loss = sapply(tree$info$splitpath, function(x) x$loss),
             npar = sapply(tree$info$splitpath, function(x) x$npar),
             nspl = sapply(tree$info$splitpath, function(x) x$nspl))
    
    ## maximum dfsplit such that the larger model improves the fit
    tab$dfsplit <- c(-(diff(tab$loss) + dfpar * diff(tab$npar)) / diff(tab$nspl), 0)
    subs <- tab$dfsplit > dfsplit
    
    if (any(subs)) {
      
      maxstep <- max(tab$step[subs])
      tree <- prune(tree, maxstep = maxstep)
      tab <- tab[subs, ]
    }
    tree$info$prunepath <- tab
  }
  
  ## return pruned model  
  return(tree)
}

prunepath.tvcm <- function(tree, steps = 1L, ...) {
  rval <- tree$info$prunepath[steps]
  class(rval) <- "prunepath.tvcm"
  return(rval)
}

print.prunepath.tvcm <- function(x, ...) {
  for (i in seq_along(x)) {
    if (!is.null(x[[i]])) {
      if (i != 1L) cat("\n")
      cat("Step:", x[[i]]$step, "\n")
      print(x[[i]]$tab, ...)
    }
  }
}

summary.tvcm <- function(object, ...)
  tvcm_print(object, type = "summary", ...)

ranef.tvcm <- function(object, ...)
  return(ranef(object$info$model, ...))

resid.tvcm <- function(object, ...)
  return(resid(object = object$info$model, ...))

residuals.tvcm <- resid.tvcm

splitpath.tvcm <- function(tree, steps = 1L,
                           details = FALSE, ...) {
  
  rval <- tree$info$splitpath[steps]
  if (!details) {
    for (i in seq_along(steps)) {
      rval[[i]]$sctest <- NULL
      rval[[i]]$lossgrid <- NULL
    }
  }
  class(rval) <- "splitpath.tvcm"
  return(rval)
}


print.splitpath.tvcm <- function(x, ...) {  
  for (i in seq_along(x)) {
    if (i != 1L) cat("\n")
    cat("Step:", x[[i]]$step)
    if (is.null(unlist(x[[i]]$varid))) {
      cat(" (no splitting processed)\n")
    } else {
      cat("\nPartition:", LETTERS[x[[i]]$partid])
      cat("\nSplitting variable:", x[[i]]$var)
      cat("\nNode:", x[[i]]$node)
      cat("\nCutpoint: ")
      cat(paste("{",paste(x[[i]]$cutpoint, collapse = "}, {"), "}\n", sep = ""))
  }
    
    if (!is.null(x[[i]]$sctest) | !is.null(x[[i]]$lossgrid))
      cat("\nDetails:\n")
    if (!is.null(x[[i]]$sctest)) {
      cat("\nCoefficient constancy tests (p-value):\n")
      for (pid in seq_along(x[[i]]$sctest)) {
        cat(paste("\nPartition ", LETTERS[pid], ":\n", sep = ""))
        print(as.data.frame(x[[i]]$sctest[[pid]]), ...)
      }
    }
    if (!is.null(x[[i]]$lossgrid)) {
      cat("\nLoss-minimizing grid search (loss reduction):\n")
      for (pid in seq_along(x[[i]]$lossgrid)) {
        for (nid in seq_along(x[[i]]$lossgrid[[pid]])) {
          for (vid in seq_along(x[[i]]$lossgrid[[pid]][[nid]])) {
            if (length(x[[i]]$lossgrid[[pid]][[nid]][[vid]]) > 0L) {
              cat("\nPartition:", names(x[[i]]$lossgrid)[pid],
                  "Node:", sub("Node", "",names(x[[i]]$lossgrid[[pid]])[nid]),
                  "Variable:", names(x[[i]]$lossgrid[[pid]][[nid]])[vid], "\n")
              print(as.data.frame(x[[i]]$lossgrid[[pid]][[nid]][[vid]]), ...)
            }
          }
        }
      }
    }
  }
}

weights.tvcm <- function(object, ...) {
  weights(extract(object, "model"))
}
