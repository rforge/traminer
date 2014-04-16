## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2014-01-03
##
## Description:
## Workhorse functions for the 'tvcolmm' function
##
## Overview:
##
## tvcolmm_fit_model:       Fit the current model
## tvcolmm_refit_model:     Refit model
## tvcolmm_fit_sctest:    Run coefficient constancy tests
## tvcolmm_fit_splitnode:   Split in variable x.
## tvcolmm_formula:         Extract separate formulas for
##                          model and partitioning from
##                          input formula.
## tvcolmm_modify_catpreds: Modify labels of categorical predictors
## tvcolmm_modify_control:  Modify the control argument before
##                          fitting the tree.
## tvcolmm_modify_modargs:  Modify the list of model args that
##                          is used for updating the model
## tvcolmm_get_start:       Get start values during partitioning
## tvcolmm_get_part:        Extract selected variables
## tvcolmm_get_terms:       Get names of varying coefficients
## tvcolmm_fit_getlevels:
## tvcolmm_fit_checksplit:  Checks whether a split is valid with
##                          respect to the minsplit parameter.
## tvcolmm_fit_surrogate:   Extracts surrogate splits using the
##                          .csurr function of the partykit
##                          package
## tvcolmm_fit_probsplit:   Extracts the conditional proportions
##                          for daughter nodes for probability
##                          splits of missing values in the
##                          splitting covariate
## tvcolmm_prune_node:
## tvcolmm_prune_depth:
## tvcolmm_prune_minsplit:
## tvcolmm_prune_alpha:
## tvcolmm_prune_nselect:
##
## Last modifications:
## 2014-04-01:   rename 'fluctest' to 'sctest'
## 2013-12-02:   remove 'tvcolmm_fit_setupnode'
## 2013-11-01:   modify 'restricted' and 'terms' correctly in
##               'tvcolmm_modify_modargs'
##
## To do:
## - improve minsplit pruning
## --------------------------------------------------------- #

tvcolmm_fit_model <- function(formula, args, control, verbose = FALSE) {

  ## set formula
  args$formula <-
    if (nlevels(args$data$Part) == 1L) formula$root else formula$tree
  
  ## fit model
  object <- try(do.call("olmm", args), TRUE)
  
  if (verbose) {
    if (nlevels(args$data$Part) == 1L) {
      cat("\n\n")
      print(object)
    } else {
      cat("\n\nPartition coefficient(s) of global node model:\n")
      terms <- if (nlevels(args$data$Part) == 1L) {
        control$terms$root
      } else {
        grep("Part", names(coef(object)), value = TRUE)
      }
      print(data.frame(Estimate = coef(object)[terms]), digits = 2)
    }
  }
  
  ## return model
  return(object)
}

tvcolmm_refit_model <- function(object, args) {

  ## set new partition
  object@frame$Part <- args$data$Part

  form <- olmm_formula(object@formula)
  contrasts <- object@contrasts

  ## update model matrix
  object@X <- olmm_mergeMm(x = model.matrix(terms(terms(form$fixefEtaVar, keep.order = TRUE)), object@frame, contrasts[intersect(names(contrasts), all.vars(form$fixefEtaVar))]), y = model.matrix(terms(form$fixefEtaInv, keep.order = TRUE), object@frame, contrasts[intersect(names(contrasts), all.vars(form$fixefEtaInv))]), TRUE)

  ## reset scores
  object@score_sbj[] <- object@score_obs[] <- object@score[] <-
    object@logLik_sbj[] <- object@logLik[] <- 0

  ## prepare optimization
  optim <- object@optim
  optim[[1]] <- object@coefficients
  optim[[4L]] <- object@restricted
  environment(optim[[2]]) <- environment()
  if (!object@dims["numGrad"]) environment(optim[[3]]) <- environment()
  if (optim$method %in% c("nlminb", "ucminf")) {
      FUN <- optim$method
      optim <- optim[-which(names(optim) == "method")] 
  } else {
      FUN <- "optim"
  }

  ## run optimization
  object@output <- do.call(FUN, optim)
  .Call("olmm_update_marg", object, object@output$par, PACKAGE = "vcolmm")

  ## return model
  return(object)
}

tvcolmm_fit_sctest <- function(model, nodes, partvar, control,
                               formula, args) {

  subject <- model@subject
  Part <- factor(fitted_node(nodes, partvar))  
  w <- weights(model)
  
  ## get column names of scores to test
  terms <- if (depth(nodes) == 0L) control$terms$root else control$terms$tree

  ## get variable types
  functional <- sapply(partvar, function(x) {
    switch(class(x)[1],
           factor = control$functional.factor,
           ordered = control$functional.ordered,
           integer = control$functional.numeric,
           numeric = control$functional.numeric)
  })
  
  
  ## prepare list with arguments for 'sctest'
  dn <- list(paste("Part", levels(Part), sep = ""), colnames(partvar))
  rval <- list(statistic = matrix(, nlevels(Part), ncol(partvar),
                 dimnames = dn),
               p.value = matrix(, nlevels(Part), ncol(partvar),
                 dimnames = dn),
               method = "M-fluctuation test")

  ## get depths of nodes for the check
  depth <- unlist(nodeapply(nodes, nodeids(nodes, terminal = TRUE),
                            function(node) info_node(node)$depth))

  ## extract transformed scores
  estfunArgs <- control$estfun
  estfunArgs$x <- model
  scores <- do.call("estfun.olmm", estfunArgs)
  
  ## hack^1: if intercept = "po" fit a second model with the second level of
  ## 'Part' as reference level
  if (nlevels(args$data$Part) > 1L && control$intercept == "po") {
    args$contrasts$Part <- contr.treatment(levels(args$data$Part), base = 2)  
    mTmp <- tvcolmm_fit_model(formula, args, control, FALSE)
    estfunArgs$x <- mTmp
    sTmp <- do.call("estfun.olmm", estfunArgs)
    subs1 <- colnames(scores) == paste("Part", levels(args$data$Part)[2], sep = "")
    subs2 <- colnames(sTmp) == paste("Part", levels(args$data$Part)[1], sep = "")
    attr <- attributes(scores)
    scores <- cbind(scores[, 1:(which(subs1) - 1) ,drop = FALSE],
                    sTmp[, which(subs2) ,drop = FALSE],
                    scores[, which(subs1):ncol(scores) ,drop = FALSE])
    attributes(scores) <- appendDefArgs(attributes(scores), attr)
    rm(mTmp, sTmp)
  }

  gefpArgs <- appendDefArgs(control$gefp, list(center = TRUE, silent = TRUE))
  gefpArgs$object <- model; gefpArgs$scores <- scores;
  
  ## apply test for each variable and partition separately  
  for (i in 1:ncol(partvar)) {    
    for (j in 1:nlevels(Part)) {
      
      ## extract observations of current partition
      rows <- Part == levels(Part)[j]
        
      ## extract variable to test
      z <- partvar[, i]
      
      ## check if partition contains at least one permitted split
      ## this ensures that we don't select a variable without permitted split
      if (is.numeric(z)) {
        sz <- cumsum(tapply(w[rows], z[rows], sum))
        mb <- max(c(0, apply(cbind(sz, sum(w[rows]) - sz), 1, min)))
      } else {
        sz <- tvcolmm_fit_getlevels(z, rows)
        mb <- max(c(0, apply(sz, 1, function(x) {
          subs <- z %in% levels(z)[x]
          min(sum(w[subs & rows]), sum(w[!subs & rows]))
        })))
      }
      if (sum(w[rows]) >= control$minsplit & 
          mb >= control$minbucket &
          depth[j] < control$maxdepth) {
      
        ## function to extract scores
        cols <- sub("Node", levels(Part)[j], terms, fixed = TRUE)
  
        ## run the tests
        gefpArgs$order.by <- z; gefpArgs$terms <- cols; gefpArgs$subset <- rows
        gefp <- try(do.call("gefp.olmm", gefpArgs), TRUE)
                                    
        if (!inherits(gefp, "try-error")) {

          ## parameters for categorical variables
          trim <- control$minbucket / sum(w[rows])
          order.by <- z[rows]
          if (is.character(functional)) {
            functional <- tolower(functional)
            fi <- switch(functional[i],
                         "suplm" = supLM(from = trim),
                         "lmuo" = catL2BB(gefp),
                         stop("Unknown efp functional."))
          } else {
            fi <- functional[i]
          }
          test <- try(sctest(gefp, functional = fi), silent = TRUE)
        } else {
          test <- gefp
        }
        
        if (!inherits(test, "try-error")) {
          
          ## extract information from test
          rval$statistic[j, i] <- test$statistic
          rval$p.value[j, i] <- test$p.value
        }
      }
    }
  }

  if (control$verbose && all(depth) == control$maxdepth)
    cat("\n\nMaximal depth reached. Return object.")
  
  ## select variables randomly (for random forests)
  valid.pval <- apply(rval$p.value, 2, function(x) sum(!is.na(x)) > 0)
  if (control$mtry < ncol(partvar) & sum(valid.pval) > 1L) {
    mtry <- min(sum(valid.pval), if (control$mtry < 0) { sample(1:sum(valid.pval), 1) } else if (control$mtry < 1) { ceiling(sum(valid.pval) * control$mtry) } else { control$mtry})
    rselect <- sort(sample(which(valid.pval), mtry))
    rval$p.value[, setdiff(1:ncol(rval$p.value), rselect)] <- NA
  }

  ## apply partial Bonferroni correction
  if (control$bonferroni) {
    for (i in 1:nlevels(Part)) {
      pval1 <- pmin(1, sum(!is.na(rval$p.value[i, ])) * rval$p.value[i, ])
      pval2 <- 1 - (1 - rval$p.value[i, ])^sum(!is.na(rval$p.value[i, ]))
      rval$p.value[i, ] <- ifelse(!is.na(rval$p.value[i,]) &
                                  (rval$p.value[i, ] > 0.01), pval2, pval1)
    }
  }

  ## print results
  if (control$verbose) {

    ## results
    cat("\nCoefficient constancy tests (p-value):\n")
    print(data.frame(rbind(
                       functional = functional,
                       format(rval$p.value, digits = 2))))

    ## add notes if necessary
    if (attr(scores, "conv.T") == 0L)
      cat("\nNote: Computation of transformation matrix failed.",
          "\nTests are based on raw scores and may be inaccurate.\n")

    if (attr(scores, "conv.complete") == 0L)
      cat("\nNote: Imputation method for complete data failed.",
          "\nThe p-values may be different in a replication.\n")
    
  }
  return(rval)
}


tvcolmm_fit_loss <- function(varid, partid, partvar, nodes,
                             model, formula, args, control, step) {

  if (control$verbose) cat("\nPerform split? ")

  ## if no splitting is necessary ...
  if (length(partid) == 0L) {
    if (control$verbose) cat("No\nToo few observations in nodes or/and",
                             "maximum depth reached. Return object.\n")
    return(NULL)
  }

  ## ... otherwise start splitting
  if (control$verbose) cat("Yes\n\n* computing splits ")
  
  Part <- factor(fitted_node(nodes, partvar))    
  w <- weights(model)
  
  ## function to get all cutpoints of a variable in a partition
  getCutpoints <- function(vid, pid) {
    subscripts <- Part == levels(Part)[pid]
    if (sum(subscripts) < 1L | sum(w[subscripts]) < control$minsplit)
      return(matrix(,0,1))
    z <- partvar[, vid]
    if (is.numeric(z)) {
      sz <- unique(sort(z[subscripts]))
      if ((length(sz) - 1)  > control$maxevalsplit) {
        rval <-  as.double(quantile(sz, ppoints(control$maxevalsplit)))
      } else {
        rval <- sz[-length(sz)]
      }
      rval <- matrix(rval, ncol = 1,
                     dimnames = list(1:length(rval), "cut"))
    } else {
      rval <- tvcolmm_fit_getlevels(z, subscripts)
    }
    return(rval)
  }

  ## function to get the splitstatistic for a given cutpoint
  getSplitstat <- function(cutpoint, vid, pid, args, model) {
    if (!(vid %in% varid & pid %in% partid)) return(Inf)
    subscripts <- Part == levels(Part)[pid]
    z <- partvar[, vid]    
    if (is.numeric(z)) {
      zs <- z <= cutpoint
    } else {
      zs <- z %in% levels(z)[cutpoint]            
    }
    args$data$Part[subscripts] <- levels(Part)[pid]
    args$data$Part[zs & subscripts] <- "Right"
    if (sum(w[zs & subscripts]) < control$minbucket) return(Inf)
    if (sum(w[!zs & subscripts]) < control$minbucket) return(Inf)
    model <- tvcolmm_refit_model(model, args)
    rval <- ifelse(inherits(model, "try-error"), NaN, control$lossfun(model))
    if (control$verbose)
      cat(if (is.nan(rval)) "f" else if (rval == Inf) "i" else ".")
    return(rval)
  }

  ## setup a prototyp model to use "tvcolmm_refit_model" function
  args$data$Part <- Part
  levels(args$data$Part) <- c(levels(args$data$Part), "Right")
  subscripts <- args$data$Part == levels(args$data$Part)[partid[1]]
  args$data$Part[sample(which(subscripts), as.integer(sum(subscripts)/2))] <- "Right"
  args$doFit <- FALSE
  splitModel <- tvcolmm_fit_model(formula, args, control)

  ## set start values
  if (control$fast > 0L) {
    subs <- intersect(names(coef(splitModel)), names(coef(model)))
    splitModel@coefficients[subs] <- model@coefficients[subs]
  }
    
  ## run computation (slow!!!)
  lossgrid <- lapply(1:ncol(partvar), function(vid) {
    lapply(1:nlevels(Part), function(pid, vid) {
      cp <- getCutpoints(vid, pid)
      st <- apply(cp, 1, getSplitstat, vid = vid, pid = pid,
                  args = args, model = splitModel)
      return(cbind(cp, loss = st))
    }, vid = vid)
  })

  Part <- factor(fitted_node(nodes, partvar))
  loss <- lapply(lossgrid, lapply, function(x) x[, ncol(x)])
  
  ## get split
  minLoss <- function(x) {
    x <- unlist(x)
    if (length(x) > 0) min(x, na.rm = TRUE) else Inf
  }
  vid <- (1:ncol(partvar))[which.min(sapply(loss, minLoss))]
  pid <- (1:nlevels(Part))[which.min(sapply(loss[[vid]], minLoss))]
  stat <- lossgrid[[vid]][[pid]]
  cid <- which.min(stat[, ncol(stat)])
  
  if (control$verbose) cat(" OK")
  
  return(list(varid = vid, partid = pid, cutid = cid, lossgrid = lossgrid))
}

tvcolmm_fit_splitnode <- function(nodes, model, loss, partvar, step) {

  Part <- factor(fitted_node(nodes, partvar))
  vid <- loss$varid
  pid <- loss$partid
  pidLab <- nodeids(as.partynode(nodes), terminal = TRUE)[pid]
  stat <- loss$lossgrid[[vid]][[pid]]
  stat <- stat[loss$cutid, ]
  cut <- stat[1:(length(stat) - 1)]
  x <- partvar[, vid]
  
  ## collect information for the split
  subscripts <- Part == levels(Part)[pid]
  if (is.numeric(x)) { # numerical variables
    breaks <- as.double(cut)
    index <- NULL
    ordered <- TRUE
    subsLeft <- subscripts & x <= breaks
    subsRight <- subscripts & x > breaks
  } else {
    subsLeft <- subscripts & x %in% levels(x)[cut == 1L]
    subsRight <- subscripts & x %in% levels(x)[cut == 0L]
    if (is.ordered(x)) { 
      breaks <- as.double(max(which(cut == 1)))
      index <- NULL
      ordered <- TRUE
    } else {
      breaks <- NULL
      index <- as.integer(-cut + 2)
      index[table(x[subscripts]) == 0] <- NA
      ordered <- FALSE
    }
  }
    
  ## get current nodes
  nodes <- as.list(nodes)
  
  ## setup 'newnodes' object
  subs <- which(sapply(nodes, function(node) node$id) == pidLab)
  newnodes <- vector("list", length(nodes) + 2)
  newnodes[1:subs] <- nodes[1:subs]
  if (length(nodes) > pidLab)
    newnodes[(subs + 3L):length(newnodes)] <-
      nodes[(subs + 1L):length(nodes)]
  
  ## adjust ids of children
  ids <- sapply(newnodes, function(x) if (!is.null(x$id)) x$id else -1)
  for (i in 1L:length(newnodes))
    if (!is.null(newnodes[[i]]$kids))
      newnodes[[i]]$kids <- which(ids %in% newnodes[[i]]$kids)
  
  ## setup new split
  newnodes[[subs]]$split <-
    partysplit(varid = vid, breaks = breaks, index = index,
               info = list(ordered = ordered, step = step))
  newnodes[[subs]]$kids <- pidLab + 1L:2L
  newnodes[[subs]]$info$dims <-
    c(N = length(unique(model@subject[subscripts])),
      n = sum(subscripts))
  
  ## add new children
  newnodes[[subs + 1L]] <-
    list(id = pidLab + 1L,
         info = list(
           dims = c(
             N = length(unique(model@subject[subsLeft])),
             n = sum(subsLeft)),
           depth = newnodes[[subs]]$info$depth + 1L))
  
  newnodes[[subs + 2L]] <-
    list(id = pidLab + 2L,
         info =
         list(dims = c(
                N = length(unique(model@subject[subsRight])),
                n = sum(subsRight)),
              depth = newnodes[[subs]]$info$depth + 1L))
  
  ## adjust ids
  for (i in 1L:length(newnodes))
    newnodes[[i]]$id <- i
  
  ## return new nodes
  return(as.partynode(newnodes))
}

tvcolmm_formula <- function(model, control, env = parent.frame()) {

  terms <- control$terms$root
  terms <- unique(sub("Eta[1-9]+:", "", terms))  
  intercept <- control$intercept
  
  ## modify predictor-variable terms
  fixefEtaVar <- names(fixef(model, "npo"))
  fixefEtaVar <- unique(sub("Eta([1-9]?[0-9]?[0-9]):", "", fixefEtaVar))
  fixefEtaVar <- setdiff(fixefEtaVar, "(Intercept)")
  fixefEtaVar <- gsub("[[:punct:]]", "", fixefEtaVar)
  if (any(fixefEtaVar %in% terms))
    fixefEtaVar[fixefEtaVar %in% terms] <-
      paste("Part:", fixefEtaVar[fixefEtaVar %in% terms], sep = "")
  if (intercept == "npo")
    fixefEtaVar <- c("Part", fixefEtaVar, "-1")
  if (length(fixefEtaVar) == 0) fixefEtaVar <- c("1")
  
  ## modify predictor-invariant terms
  fixefEtaInv <- names(fixef(model, "po"))
  fixefEtaInv <- gsub("[[:punct:]]", "", fixefEtaInv)
  if (any(fixefEtaInv %in% terms))
    fixefEtaInv[fixefEtaInv %in% terms] <-
      paste("Part:", fixefEtaInv[fixefEtaInv %in% terms], sep = "")
  if (intercept == "po")
    fixefEtaInv <- c("Part", fixefEtaInv)
  
  ## get response variable name
  yName <- all.vars(formula(model))[1]

  W <- model.matrix(model, "ranef")
  
  ## predictor-variable random effects
  if (model@dims["hasRanef"] > 0) {
    ranefEtaVar <- colnames(W, "ranef")[attr(W, "merge") == 1L]
    ranefEtaVar[ranefEtaVar == "(Intercept)"] <- "1"
    ranefEtaVar <- gsub("[[:punct:]]", "", ranefEtaVar)
    if (length(ranefEtaVar) > 0) {   
      ranefEtaVar <- paste("(", paste(ranefEtaVar, collapse = " + "), " | ",
                           model@subjectName, ")", sep = "")
    } 
  } else {
    ranefEtaVar <- c()
  }
  
  ## predictor-invariant random effects
  if (model@dims["hasRanef"] > 0) {
    ranefEtaInv <- colnames(W, "ranef")[attr(W, "merge") == 2L]
    ranefEtaInv[ranefEtaInv == "(Intercept)"] <- "1"
    ranefEtaInv <- gsub("[[:punct:]]", "", ranefEtaInv)
    if (length(ranefEtaInv) > 0) {   
      ranefEtaInv <- paste("(", paste(ranefEtaInv, collapse = " + "), " | ",
                           model@subjectName, ")", sep = "")
      
    }
  } else {
    ranefEtaInv <- c()
  }

  ## paste together the formula
  return(formula(paste(yName, " ~ ",
                       paste(c(fixefEtaInv, ranefEtaInv),
                               collapse = " + "), " | ",
                       paste(c(fixefEtaVar, ranefEtaVar),
                             collapse = " + "), sep = ""),
                 env = environment(formula(model))))
}


tvcolmm_modify_catpreds <- function(formula, data) {

  dataClasses <- unlist(lapply(data, function(x) class(x)[1]))
  form <- olmm_formula(formula)
  modelVars <- setdiff(all.vars(form$full)[-1], form$subjectName)
  for (i in which(dataClasses %in% c("ordered", "factor") &
                  colnames(data) %in% modelVars)) {
    levs <- levels(data[, i])
    levels(data[, i]) <- gsub(" +", "", levels(data[, i]))
    levels(data[, i]) <- gsub("[[:punct:]]", "", levels(data[, i]))
    if (!identical(levs, levels(data[, i])))
      warning(paste("levels of variable '",
                    colnames(data)[i],
                    "' were modified.", sep = ""))
    
  }
  return(data)
}

  
tvcolmm_modify_control <- function(model, control) {

  ## prevent conflicts
  if (control$intercept == "none") {
    control$restricted <- unique(c("(Intercept)", control$restricted))
    control$terms <- setdiff(control$terms, "(Intercept)")
    intercept <- FALSE
  } else {
    intercept <- TRUE
  }
  
  ## get needed information
  atEtaVar <- terms(model, "fixef-npo")
  atEtaInv <- terms(model, "fixef-po")
  atEtaVar <- c("(Intercept)", attr(atEtaVar, "term.labels"))
  atEtaInv <- attr(atEtaInv, "term.labels")
  merge <- attr(model.matrix(model), "merge")
  assignEtaVar <- attr(model.matrix(model), "assign")[merge == 1L]
  assignEtaInv <- attr(model.matrix(model), "assign")[merge == 2L]
  merge <- attr(model.matrix(model), "merge")
  fixefEtaVar <- names(fixef(model, "npo"))
  fixefEtaInv <- names(fixef(model, "po"))
  mmNamesEtaVar <- colnames(model.matrix(model))[merge == 1L]
  mmNamesEtaInv <- colnames(model.matrix(model))[merge == 2L]

  ## objects to determine
  termsFixefEtaVar <- termsFixefEtaInv <-
    restFixefEtaVar <- restFixefEtaInv <- character()
  
  ## modify entries in 'terms'
  if (length(control$terms) > 0L) {

    if (all(control$terms %in% c(fixefEtaVar, fixefEtaInv))) {

      termsFixefEtaVar <- intersect(control$terms, fixefEtaVar)
      termsFixefEtaVar <- unique(sub("Eta[1-9]+:", "", termsFixefEtaVar))
      termsFixefEtaInv <- intersect(control$terms, fixefEtaInv)
      
    } else if (all(control$terms %in% c(atEtaVar, atEtaInv))) {

      subs <- which(atEtaVar %in% control$terms)
      if (length(subs) > 0)
        termsFixefEtaVar <- mmNamesEtaVar[assignEtaVar %in% (subs - 1L)]
      
      subs <- which(atEtaInv %in% control$terms)
      if (length(subs) > 0)
        termsFixefEtaInv <- mmNamesEtaInv[assignEtaInv %in% subs]
      
    } else {

      stop("invalid 'terms' argument")
    }
    if (!"(Intercept)" %in% c(termsFixefEtaVar, termsFixefEtaVar))
      intercept <- FALSE
  }

  if (length(control$restricted) > 0L) {

    if (all(control$restricted %in% c(fixefEtaVar, fixefEtaInv))) {

      restFixefEtaVar <- intersect(control$restricted, fixefEtaVar)
      restFixefEtaVar <- unique(sub("Eta[1-9]+:", "", restFixefEtaVar))
      restFixefEtaInv <- intersect(control$restricted, fixefEtaInv)
      
    } else if (all(control$restricted %in%
                   c(atEtaVar, atEtaInv))) {

      subs <- which(atEtaVar %in% control$restricted)
      if (length(subs) > 0)
        restFixefEtaVar <- mmNamesEtaVar[assignEtaVar %in% (subs - 1)]
      
      subs <- which(atEtaInv %in% control$restricted)
      if (length(subs) > 0)
        restFixefEtaInv <- mmNamesEtaInv[assignEtaInv %in% subs]
      
    } else {

      stop("invalid 'restricted' argument")
    }
  } 

  ## if 'terms' is empty, add all fixed effect parameters except
  ## 'restricted' parameters
  if (length(control$terms) == 0L) {

    termsFixefEtaVar <-
      setdiff(unique(sub("Eta[1-9]+:", "", fixefEtaVar)), restFixefEtaVar)
    termsFixefEtaInv <- setdiff(fixefEtaInv, restFixefEtaInv)
  }

  if (!control$intercept == "npo")
    termsFixefEtaVar <- setdiff(termsFixefEtaVar, "(Intercept)")

  ## if 'restricted' is empty, add all fixef effect parameters
  ## except 'terms' parameters
  if (length(control$restricted) == 0L) {

    restFixefEtaVar <-
      setdiff(unique(sub("Eta[1-9]+:", "", fixefEtaVar)), termsFixefEtaVar)
    restFixefEtaInv <-
      setdiff(fixefEtaInv, termsFixefEtaInv)
  }

  ## erase ':' and several other signs
  FUN1 <- function(x) {
    rval <- x
    subsInt <- rval == "(Intercept)"
    rval[!subsInt] <- gsub("[[:punct:]]", "", x[!subsInt])
    if (any(!x %in% rval)) {
      subs <- !x %in% rval
      warning(paste("term(s)", paste("'", x[subs],"'", collapse = ",", sep = ""), "were changed to", paste("'", rval[subs],"'", collapse = ",", sep = "")))
    }
    return(rval)
  }
  
  termsFixefEtaVar <- FUN1(termsFixefEtaVar)    
  termsFixefEtaInv <- FUN1(termsFixefEtaInv)
  restFixefEtaVar <- FUN1(restFixefEtaVar)
  restFixefEtaInv <- FUN1(restFixefEtaInv)
  
  ## append 'Eta[1-9]+':
  FUN2 <- function(x, n)
    paste("Eta", rep(1:n, each = length(x)), ":", rep(x, n), sep = "")
  if (length(termsFixefEtaVar) > 0)
    termsFixefEtaVar <-
      FUN2(termsFixefEtaVar, model@dims["nEta"])
  if (length(restFixefEtaVar) > 0)
    restFixefEtaVar <-
      FUN2(restFixefEtaVar, model@dims["nEta"])

  control$terms <- list(original = c(unique(sub("Eta[1-9]+:", "", termsFixefEtaVar)), termsFixefEtaInv), root = c(termsFixefEtaVar, termsFixefEtaInv))
  
  FUN3 <- function(x) {
    if (grepl("(Intercept)", x)) {
      rval <- gsub("(Intercept)", "PartNode", x, fixed = TRUE) 
    } else {
      if (grepl("Eta[1-9]+:", x)) {
        rval <- strsplit(x, ":")[[1]]
        rval <- paste(rval[1], "PartNode", rval[-1], sep = ":")
      } else {
        rval <- paste("PartNode:", x, sep = "")
      }
    }
    return(rval)
  }
  termsFixefEtaVar <- unlist(lapply(termsFixefEtaVar, FUN3))
  termsFixefEtaInv <- unlist(lapply(termsFixefEtaInv, FUN3))

  if (intercept && control$intercept == "po")
    termsFixefEtaInv <- c("PartNode", termsFixefEtaInv)
  
  control$terms$tree <- c(termsFixefEtaVar, termsFixefEtaInv)  
  control$restricted <- c(setdiff(unique(sub("Eta[1-9]+:", "", restFixefEtaVar)), "(Intercept)"), restFixefEtaInv)
  
  ## coefficients to treated as nuisance for estfun.olmm
  if (is.null(control$estfun$nuisance)) {
    nui <- control$restricted
    if (model@dims["hasRanef"] > 0)
      nui <- c(nui, grep("ranefCholFac", names(coef(model)), value = TRUE))
    control$estfun$nuisance <- nui
  }  
  
  if (length(control$terms$root) * length(control$terms$tree) == 0)
    stop("no 'terms' found.")
  
  return(control)
}


tvcolmm_modify_modargs <- function(model, args, control) {

  yName <- all.vars(formula(model))[1]
  subjectName <- model@subjectName
  
  ## modify formulas for fixed effects
  X <- model.matrix(model, "fixef")
  xTerms <- colnames(X)
  subsInt <- xTerms == "(Intercept)"
  xTerms[!subsInt] <- gsub("[[:punct:]]", "", xTerms[!subsInt])
  colnames(X) <- xTerms

  if (any(subs <- (attr(X, "merge") == 1 & attr(X, "assign") > 0))) {
    fixefEtaVar <- paste(xTerms[subs], collapse = " + ")
  } else {
    fixefEtaVar <- "1"
  }

  if (any(subs <- (attr(X, "merge") == 2 & attr(X, "assign") > 0))) {
    fixefEtaInv <- paste(xTerms[subs], collapse = " + ")
  } else {
    fixefEtaInv <- "-1"
  }

  ## modify formulas for random effects
  W <- model.matrix(model, "ranef")
  wTerms <- colnames(W)
  subsInt <- wTerms == "(Intercept)"
  wTerms[!subsInt] <- gsub("[[:punct:]]", "", wTerms[!subsInt])
  colnames(W) <- wTerms
  if (model@dims["hasRanef"] > 0) {
    if ("(Intercept)" %in% wTerms)
      wTerms[wTerms == "(Intercept)"] <- "1"
    
    if (any(subs <- (attr(W, "merge") == 1 & attr(W, "assign") >= 0))) {
      ranefEtaVar <- paste(wTerms[subs], collapse = " + ")
      ranefEtaVar <- paste("(", ranefEtaVar, "|", subjectName, ")", sep = "")
    } else {
      ranefEtaVar <- ""
    }
    
    if (any(subs <- (attr(W, "merge") == 2 & attr(W, "assign") >= 0))) {
      ranefEtaInv <- paste(wTerms[subs], collapse = " + ")
      ranefEtaInv <- paste("(", ranefEtaInv, "|", subjectName, ")", sep = "")
    } else {
      ranefEtaInv <- ""
    }
  } else {
    ranefEtaVar <- ranefEtaInv <- ""
  }

  ## paste all together to one formula
  ff <- formula(paste(yName, " ~ ",
                      fixefEtaInv, ifelse(ranefEtaInv != "",  " + ", ""),
                      ranefEtaInv, " | ",
                      fixefEtaVar, ifelse(ranefEtaVar != "",  " + ", ""),
                      ranefEtaVar, sep = ""))
  environment(ff) <- environment(args$formula)
  args$formula <- ff

  ## modify data
  
  args$data <- as.data.frame(X)
  args$data <-
    args$data[, colnames(args$data) != "(Intercept)", drop = FALSE]
  args$data[, c(yName, subjectName)] <-
    model.frame(model)[, c(yName, subjectName)]
  wVars <- setdiff(colnames(W), colnames(X))
  if (length(wVars) > 0)
    args$data[, wVars] <- W[, wVars]

  ## delete partitioned effects from initial values
  FUN <- function(x) {
    if (is.null(x)) return(x)
    names.only <- FALSE
    if (is.null(names(x))) {
      names(x) <- x
      names.only <- TRUE
    }
    subs <- grepl("ranefCholFac", names(x))
    xFE <- x[!subs]
    xRE <- x[subs]
    dims <- model@dims
    newNames <- c(paste("Eta", rep(1:dims["nEta"], each = dims["pEtaVar"]), ":", rep(xTerms[1:dims["pEtaVar"]], dims["nEta"]), sep = ""), xTerms[(dims["pEtaVar"] + 1):dims["pEta"]])
    names(newNames) <- names(fixef(model))
    xFE <- xFE[intersect(names(newNames), names(xFE))]
    names(xFE) <- newNames[intersect(names(newNames), names(xFE))]
    rval <- c(xFE, xRE)    
    if (names.only) { 
      return(names(rval))
    } else {
      return(rval)
    }
  }
  args$start <- FUN(args$start)
  args$restricted <- FUN(args$restricted)
  
  return(args)
}

tvcolmm_get_start <- function(model, args) {
  coef <- coef(model)
  FUN <- function(x) {
    if (!grepl("Part", x)) {
      return(TRUE)
    } else {
      tmp <- paste("Part", levels(args$data$Part), sep = "")
      return(any(sapply(tmp, grepl, x)))
    }
  }
  return(coef[sapply(names(coef), FUN)])
}

tvcolmm_get_part <- function(object, data, weights = rep(1.0, nrow(data))) {
  fitted <-
    fitted_node(object$node, data[, colnames(object$data), drop = FALSE])
  rval <- factor(fitted, levels = nodeids(object$node, terminal = TRUE))
  if (nlevels(rval) > 1L) {
    con <- contr.sum(levels(rval))
    tab <- tapply(weights, rval, sum)
    con[nrow(con),] <- con[nrow(con),] * tab[-length(tab)] / tab[length(tab)]
    colnames(con) <- levels(rval)[1:(nlevels(rval) - 1)]
    contrasts(rval) <- con
  }
  return(rval)
}

tvcolmm_get_terms <- function(object) {
  rval <- NULL
  ids <- nodeids(object)
  model <- extract(object, "model")
  
  formula <- olmm_formula(object$info$formula$tree)$fixefEtaVar 
  factors <- attr(terms(formula), "factors")
  if ("Part" %in% all.vars(formula)) {    
    tmp1 <- paste("Eta", 1:model@dims["nEta"], ":", sep  = "")
    tmp2 <- colnames(factors)[factors["Part", ] > 0]
    rval <- c(rval, paste(rep(tmp1, each = length(tmp2)),
                          rep(tmp2, length(tmp1)), sep = ""))
    
  }
  
  formula <- olmm_formula(object$info$formula$tree)$fixefEtaInv 
  factors <- attr(terms(formula), "factors")
  if ("Part" %in% all.vars(formula)) {
    tmp <- c(colnames(factors)[factors["Part", ] > 0])
    rval <- c(rval, c(colnames(factors)[factors["Part", ] > 0]))
  }
  return(rval)
}

tvcolmm_fit_getlevels <- function(x, subscripts) {

  ## all levels
  nl <- nlevels(x)

  ## levels in current partition
  xd <- droplevels(x[subscripts > 0])
  nld <- nlevels(xd)
  
  if (inherits(x, "ordered")) {
    indx <- diag(nl)
    indx[lower.tri(indx)] <- 1
    indx <- indx[-nl, , drop = FALSE]
    rownames(indx) <- levels(x)[-nl]
  } else {
    mi <- 2^(nld - 1) - 1
    indx <- matrix(0, nrow = mi, ncol = nl)
    for (i in seq(1, mi, length.out = mi)) { # go through all splits #
      ii <- i
      for (l in 1:nld) {
        indx[i, which(levels(x) %in% levels(xd))[l]] <- ii %% 2;
        ii <- ii %/% 2   
      }
    }
    rownames(indx) <- apply(indx, 1, function(z) paste(levels(x)[!is.na(z) & z > 0], collapse = "+"))
  }
  colnames(indx) <- as.character(levels(x))
  storage.mode(indx) <- "logical"
  return(indx)
}


## tvcolmm_fit_surrogate <- function(split, partvar, weights, model, control) {

##   if (control$condsurrogate == "subject") {

##     subject <- model@subject
##     weights <- unique(cbind(subject, weights))[, -1]
    
##     FUN <- function(x) {
##       rval <- x[1]
##       x <- na.omit(x)
##       if (length(x) == 0) rval[1] <- NA else rval <- sample(x, 1)
##       return(rval)
##     }

##     partvar <- aggregate(partvar, by = list(subject), FUN = FUN)
##     partvar <- partvar[, -1, drop = FALSE]
##   }
  
##   x <- partvar[, split$varid]
##   xna <- is.na(x)
  
##   weights[xna] <- 0L
  
##   inp <- rep(TRUE, ncol(partvar))
##   inp[split$varid] <- FALSE

##   kidsids <- kidids_split(split, data = partvar)
##   kidsids[is.na(kidsids)] <- 1L

##   if (sum(inp) > 0) {

##     rval <-
##       partykit:::.csurr(split = kidsids, data = partvar,
##                         inp = inp, weights = as.integer(weights),
##                         ctrl = list(maxsurrogate = control$maxsurrogate))
##   } else {

##     rval <- list()
##   }
  
##   return(rval)
## }


## tvcolmm_fit_probsplit <- function(split, model, partvar, weights, control) {

##   if (control$probsurrogate) {
        
##     kidsids <- kidids_split(split, data = partvar)[weights > 0]    

##     ## choose randomly a value of an individual
##     if (control$condsurrogate == "subject") {
##       FUN <- function(x) {
##         rval <- x[1]
##         x <- na.omit(x)
##         if (length(x) == 0) rval[1] <- NA else rval <- sample(x, 1)
##         return(rval)
##       }
##       kidsids <- tapply(kidsids, as.integer(model@subject[weights > 0]), FUN)
##     }

##     rval <- prop.table(table(kidsids))
##     storage.mode(rval) <- "double"
##     return(rval)
    
##   } else {

##     return(NULL)
##   }
## }


tvcolmm_prune_node <- function(object, alpha = NULL,
                               depth = NULL, width = NULL,
                               minsplit = NULL, minbucket = NULL,
                               nselect = NULL, step = NULL) {

  stopifnot(class(object)[1] %in% c("tvcolmm", "party", "partynode"))
  if ("partynode" %in% class(object)) {
    rval <- object
  } else {
    rval <- object$node
  }

  control <- object$info$control
  if (!is.null(depth) && depth(rval) > 0L)
    rval <- tvcolmm_prune_depth(rval, 0L, depth)
  if (!is.null(width) && width(rval) > 0L) 
    rval <- tvcolmm_prune_step(object$node, width)
  if (!is.null(minsplit) && depth(rval) > 0L)
    rval <- tvcolmm_prune_minsplit(rval, minsplit)
  if (!is.null(minbucket) && depth(rval) > 0L)
    rval <- tvcolmm_prune_minbucket(rval, minbucket)
  if (!is.null(nselect) && length(extract(object, "selected"))) {
    splitpath <- object$info$splitpath
    vars <- sapply(splitpath, function(x) x$varid)
    nvars <- sapply(1:length(vars), function(x) length(unique(vars[1:x])))
    step <- max(c(0, which(nvars <= nselect)))
    rval <- tvcolmm_prune_step(object$node, step)
  }
  if (!is.null(alpha) && depth(rval) > 0L) {
    splitpath <- object$info$splitpath
    ids <- nodeids(object)
    stepids <- nodeapply(object,  ids, function(node) node$split$info$step)
    stepids <- sapply(stepids, function(x) if (is.null(x)) Inf else x)
    p.value <- extract(object, "p.value")
    step <- which(!is.na(p.value) & p.value <= alpha)[1]
    if (length(step) > 0L)
      rval <- tvcolmm_prune_step(object$node, step)
  }
  if (!is.null(step) && object$info$nsteps > step) {
    rval <- tvcolmm_prune_step(object$node, step)
  }
  
  ## delete empty nodes and adjust id labeling
  rval <- as.list(rval)
  rval <- rval[sapply(rval, function(x) !is.null(x))] # delete nodes
  ids_old <- unlist(lapply(rval, function(x) x$id))
  ids_new <- 1L:length(rval)
  for (i in ids_new) {
    rval[[i]]$id <- ids_new[ids_old == rval[[i]]$id]
    for (j in 1:length(rval[[i]]$kids))
      rval[[i]]$kids[j] <- ids_new[ids_old == rval[[i]]$kids[j]]
  }
  
  rval <- as.partynode(rval)
  return(rval)
}

tvcolmm_prune_depth <- function(node, d, depth) {
  
  if (!is.terminal(node)) {
    if (d + 1 <= depth) {
      kids <- sapply(node$kids, function(kids) kids$id)
      for (i in 1:length(kids)) {
        node$kids[[i]] <-
          tvcolmm_prune_depth(node$kids[[i]], d + 1L, depth)
      }
    } else {
      node$kids <- NULL
      node$surrogates <- NULL
      node$split <- NULL
    }
  }
  return(node)
}

tvcolmm_prune_minsplit <- function(node, minsplit) {
  
  if (!is.terminal(node)) {
    n <- info_node(node)$dims["n"]
    if (n >= minsplit) {
      kids <- sapply(node$kids, function(kids) kids$id)
      for (i in 1:length(kids)) {
        node$kids[[i]] <- tvcolmm_prune_minsplit(node$kids[[i]], minsplit)
      }
    } else {
      node$kids <- NULL
      node$surrogates <- NULL
      node$split <- NULL
    }
  }
  return(node)
}

tvcolmm_prune_minbucket <- function(node, minbucket) {
  
  if (!is.terminal(node)) {
    n <- sapply(node$kids, function(x) x$info$dims["n"])
    if (min(n) >= minbucket) {
      kids <- sapply(node$kids, function(kids) kids$id)
      for (i in 1:length(kids)) {
        node$kids[[i]] <- tvcolmm_prune_minbucket(node$kids[[i]], minbucket)
      }
    } else {
      node$kids <- NULL
      node$surrogates <- NULL
      node$split <- NULL
    }
  }
  return(node)
}

tvcolmm_prune_step <- function(node, step) {
  
  if (!is.terminal(node)) {
    if (node$split$info$step <= step) {
      kids <- sapply(node$kids, function(kids) kids$id)
      for (i in 1:length(kids)) {
        node$kids[[i]] <- tvcolmm_prune_step(node$kids[[i]], step)
      }
    } else {
      node$kids <- NULL
      node$surrogates <- NULL
      node$split <- NULL
    }
  }
  return(node)
}

tvcolmm_get_estimates <- function(object, what = c("coef", "sd", "var"), ...) {

  what <- match.arg(what)
  
  ids <- nodeids(object, terminal = TRUE)
  control <- extract(object, "control")

  rval <- list()
  
  ## extract coefficients
  coef <- switch(what,
                 coef = coef(extract(object, "model")),
                 sd = diag(vcov(extract(object, "model"))),
                 var = diag(vcov(extract(object, "model"))))

  ## restricted coefficients
  rval$restricted <- coef[!grepl("Part", names(coef))]
  
  ## varying coefficients
  if (depth(object) > 0L) {
    
    ## create object with all possible coefficients for each node
    terms <- tvcolmm_get_terms(object)
    varcoef <- rep(NA, length(terms) * length(ids))
    for (i in 1:length(terms))
      for (j in 1:length(ids))
        names(varcoef)[length(ids) * (i - 1) + j] <-
          sub("Part", paste("Part", ids[j], sep = ""), terms[i])
    subs <- intersect(names(varcoef), names(coef))
    varcoef[subs] <- coef[subs]
    
    if ((subs <- paste("Part", max(ids), sep = "")) %in%
        names(varcoef)) {
      con <- extract(object, "model")@contrasts$Part

      if (what == "coef") {
        varcoef[subs] <-
          sum(con[as.character(max(ids)),] *
              coef[paste("Part", setdiff(ids, max(ids)), sep = "")])
        
      } else if (what %in% c("sd", "var")) {
        varcoef[subs] <-
          sum((con[as.character(max(ids)),])^2 *
              coef[paste("Part", setdiff(ids, max(ids)), sep = "")])
      }
    } 
    
    ## create a matrix of coefficients
    FUN <- function(i) {
      name <- paste("Part", i, sep = "")
      parts <- strsplit(names(varcoef), ":")
      subs <- sapply(parts, function(x) sum(x == name) > 0)
      rval <- varcoef[subs]
      names(rval) <- sub(name, "Part", names(rval))
      names(rval) <- sub("Part:", "", names(rval))
      return(rval)
    }
    
    rval$varying <- lapply(as.character(ids), FUN)
    rval$varying <-
      matrix(unlist(rval$varying), nrow = length(ids), byrow = TRUE,
             dimnames = list(ids, names(rval$varying[[1]])))
    
  } else {
    terms <- setdiff(tvcolmm_get_terms(object), "Part")
    terms <- sapply(terms, function(x) sub("Part:", "", x))
    if (object$info$vi == "npo") {
      nEta <- object$info$model@dims["nEta"]
      subs <- names(coef) %in% paste("Eta", 1:nEta, ":(Intercept)", sep = "")
      names(coef)[subs] <- paste("Eta", 1:nEta, ":Part", sep = "")
    }
    rval$varying <- matrix(coef[terms], 1, dimnames = list("1", terms))
    rval$restricted <- coef[setdiff(names(coef), terms)]
  }

  if (what == "sd") {
    rval$restricted <- sqrt(rval$restricted)
    rval$varying <- sqrt(rval$varying)
  }
  
  return(rval)
}

tvcolmm_modify_splitpath <- function(splitpath, nodes, partvar, control) {

  ids <- nodeids(nodes)
  stepids <- nodeapply(nodes,  ids, function(node) node$split$info$step)
  stepids <- sapply(stepids, function(x) if (is.null(x)) Inf else x)
  
  for(step in 1:length(splitpath)) {
    if (step > 1L) {
      parentids <- ids[stepids < step]
      kids <- nodeapply(nodes, parentids, function(node) node$kids)
      kidids <- c()
      for (i in 1:length(kids)) {
        kidids <- c(kidids, unlist(lapply(kids[[i]], function(x) x$id)))
      }
      kidids <- setdiff(sort(kidids), parentids)
    } else {
      kidids <- 1L
    }
    splitpath[[step]]$ids <- kidids
    if (control$sctest) {
      if (!is.null(splitpath[[step]]$sctest)) {
        rownames(splitpath[[step]]$sctest$statistic) <-
          paste("Part", kidids, sep = "")
        rownames(splitpath[[step]]$sctest$p.value) <-
          paste("Part", kidids, sep = "")
      }
      if (!is.null(splitpath[[step]]$varid)) {
        if (length(splitpath[[step]]$lossgrid) > 1L)
          splitpath[[step]]$lossgrid[-splitpath[[step]]$varid] <- NULL
        if (length(splitpath[[step]]$lossgrid[[1]]) > 1L)
          splitpath[[step]]$lossgrid[[1]][-splitpath[[step]]$partid] <- NULL
        names(splitpath[[step]]$lossgrid)[1L] <- 
          colnames(partvar)[splitpath[[step]]$varid]
        names(splitpath[[step]]$lossgrid[[1]])[1L] <- 
          paste("Part", ids[stepids == step], sep = "")
      }
    } else {
      if (!is.null(splitpath[[step]]$lossgrid)) {
        names(splitpath[[step]]$lossgrid) <- colnames(partvar)
        for (i in 1:ncol(partvar)) 
          names(splitpath[[step]]$lossgrid[[i]]) <- 
            paste("Part", kidids, sep = "")
      }
    }
    if (!is.null(splitpath[[step]]$varid)) {
      splitpath[[step]]$partid <- ids[stepids == step]
      splitpath[[step]]$var <- colnames(partvar)[splitpath[[step]]$varid]
      splitpath[[step]]$cutpoint <- character_split(nodeapply(nodes, ids[stepids == step], function(node) node$split)[[1]], data = partvar)$levels
    }
  }
  class(splitpath) <- "splitpath.tvcolmm"
  return(splitpath)
}
