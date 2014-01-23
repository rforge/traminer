## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2013-01-16
##
## Description:
## Workhorse functions for the 'tvcolmm' function
##
## Overview:
##
## tvcolmm_fit_model:       Fit the current model
## tvcolmm_fit_fluctest:    Run coefficient constancy tests
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
##                        respect to the minsplit parameter.
## tvcolmm_fit_surrogate:   Extracts surrogate splits using the
##                        .csurr function of the partykit
##                        package
## tvcolmm_fit_probsplit:   Extracts the conditional proportions
##                        for daughter nodes for probability
##                        splits of missing values in the
##                        splitting covariate
## tvcolmm_prune_node:
## tvcolmm_prune_depth:
## tvcolmm_prune_minsplit:
## tvcolmm_prune_alpha:
## tvcolmm_prune_nselect:
##
## Last modifications:
## 2013-12-02:   remove 'tvcolmm_fit_setupnode'
## 2013-11-01:   modify 'restricted' and 'terms' correctly in
##               'tvcolmm_modify_modargs'
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

  object@frame$Part <- args$data$Part
  form <- olmm_formula(object@formula)
  contrasts <- object@contrasts
  object@X <- olmm_mergeMm(x = model.matrix(terms(terms(form$fixefEtaVar, keep.order = TRUE)), object@frame, contrasts[intersect(names(contrasts), all.vars(form$fixefEtaVar))]), y = model.matrix(terms(form$fixefEtaInv, keep.order = TRUE), object@frame, contrasts[intersect(names(contrasts), all.vars(form$fixefEtaInv))]), TRUE)
    
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
  object@output <- do.call(FUN, optim)
  .Call("olmm_update_marg", object, object@output$par, PACKAGE = "vcolmm")

  ## return model
  return(object)
}

tvcolmm_fit_fluctest <- function(model, nodes, partvar, control,
                                 formula, args) {

  subject <- model@subject
  Part <- model@frame$Part
  if (is.null(Part)) Part <- factor(rep(1, nrow(partvar)))

  ## get column names of scores to test
  terms <- if (depth(nodes) == 0L) control$terms$root else control$terms$tree
    
  ## get variable types
  level <- sapply(control$type.vars,
                  function(x) ifelse(x == "tm-var", "observation", "subject"))
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
               method = "M-fluctuation test",
               data.name = "model")

  ## get depths of nodes for the check
  depth <- unlist(nodeapply(nodes, nodeids(nodes, terminal = TRUE),
                            function(node) info_node(node)$depth))

  ## hack^1: if intercept = "po" fit a second model with the second level of
  ## 'Part' as reference level
  if (nlevels(args$data$Part) > 1L && control$intercept == "po") {
    args$contrasts$Part <- contr.treatment(levels(args$data$Part), base = 2)
    model2 <- tvcolmm_fit_model(formula, args, control, FALSE)
  }
  
  ## apply test for each variable and partition separately  
  for (i in 1:ncol(partvar)) {    
    for (j in 1:nlevels(Part)) {

      ## hack^1: choose the model (only important if 'intercept == "po"')
      m <- if (nlevels(Part) > 1L && control$intercept == "po" && j == 1L) model2 else model
      
      ## extract observations of current partition
      rows <- Part == levels(Part)[j]

      ## extract variable to test
      z <- partvar[rows, i]
      if (is.factor(z)) z <- droplevels(z)
      
      ## check if partition contains permitted split
      if (sum(cumsum(sort(table(z))) >= control$minsplit) > 1 &
          max(cumsum(sort(table(z)))) >= 2 * control$minsplit &
          depth[j] < control$maxdepth) {

        ## function to extract scores
        cols <- sub("Node", levels(Part)[j], terms, fixed = TRUE)
        
        ## break the ties
        if (control$breakties) {
          oi <- sample(1:length(z), length(z))
        } else {
          oi <- 1:length(z)
        }
        
        z <- z[oi]

        ## define the score function extractor for 'sctest'
        scores <- function(x) {
          rval <- estfun(m)[rows, cols, drop = FALSE]
          return(rval[oi, , drop = FALSE])
        }

        ## covariance matrix, must be improved!!!
        vcov <- function(x, order.by, data) {
          scores <- scores(x)
          n <- nrow(scores)
          ## independence case
          v1 <- crossprod(scores / sqrt(n))
          scores <-  apply(scores, 2, tapply, droplevels(subject[rows]), sum)
          ## time-invariant
          v2 <- crossprod(scores / sqrt(n))
          rval <- (level[i] == "observation") * v1 + (level[i] == "subject") * v2
          return(rval)
        }
        
        ## arguments for test
        argNames <- names(formals(strucchange:::sctest.default))
        testArgs <- control[names(control) %in% argNames]
        testArgs <- appendDefArgs(list(x = m,
                                       order.by = z,
                                       scores = scores,
                                       sandwich = FALSE,
                                       vcov = vcov,
                                       functional = functional[i]),
                                  testArgs)
        
        ## call test
        test <- try(do.call("sctest", testArgs), silent = TRUE)
        
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
  
  ## select variables randomly
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
    cat("\nFluctuation tests (p-value):\n")
    print(data.frame(rbind(
                       type =  control$type.vars,
                       functional = functional,
                       format(rval$p.value, digits = 2))))
  }
  return(rval)
}


tvcolmm_fit_splitnode <- function(varid = 1:ncol(partvar),
                                  partid = 1:width(nodes), 
                                  partvar, nodes,
                                  model, formula, args,
                                  test, control, step) {

  if (control$verbose) cat("\nPerform split? ")

  ## if no splitting is necessary ...
  if (length(partid) == 0L) {
    if (control$verbose) cat("No\nToo few observations in nodes or/and",
                             "maximum depth reached. Return object.\n")
    return(FALSE)
  }

  ## ... otherwise start splitting
  if (control$verbose) cat("Yes\n\n* computing splits ")
  
  Part <- factor(fitted_node(nodes, partvar))    

  ## function to get all cutpoints of a variable in a partition
  getCutpoints <- function(vid, pid) {
    subscripts <- Part == levels(Part)[pid]
    x <- partvar[, vid]
    if (is.numeric(x)) {
      sx <- sort(x[subscripts])              
      sx <- sx[control$minsplit:(length(sx) - control$minsplit)]
      sx <- unique(sx)
      if ((length(sx) - 1)  > control$maxevalsplit) {
        rval <-  as.double(quantile(sx, ppoints(control$maxevalsplit)))
      } else {
        rval <- sx[-length(sx)]
      }
      rval <- matrix(rval, ncol = 1)
    } else {
      rval <- vcolmm:::tvcolmm_fit_getlevels(x, subscripts)
    }
    return(rval)
  }

  ## function to get the splitstatistic for a given cutpoint
  getSplitstat <- function(cutpoint, vid, pid, args, model) {
    if (!(vid %in% varid & pid %in% partid)) return(Inf)
    subscripts <- Part == levels(Part)[pid]
    x <- partvar[, vid]    
    if (is.numeric(x)) {
      xs <- x <= cutpoint
    } else {
      xs <- x %in% levels(x)[cutpoint]            
    }
    args$data$Part[subscripts] <- levels(Part)[pid]
    args$data$Part[xs & subscripts] <- "Right"
    if (min(table(args$data$Part)[c(levels(args$data$Part)[pid], "Right")]) <
        control$minsplit)
      return(Inf)
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
  
  ## run computation (slow!!!)
  splits <- lapply(1:ncol(partvar), function(vid) {
    lapply(1:nlevels(Part), function(pid, vid) {
      cp <- getCutpoints(vid, pid)
      st <- apply(cp, 1, getSplitstat, vid = vid, pid = pid,
                  args = args, model = splitModel)
      return(cbind(cp, stat = st))
    }, vid = vid)
  })

  loss <- lapply(splits, lapply, function(x) x[, ncol(x)])
  
  if (!any(na.omit(unlist(loss)) < Inf)) {
    if (control$verbose)
      cat("\nNo admissible split found. Return object.\n")
    return(FALSE)
  }
      
  if (control$verbose) cat(" OK")
  
  ## get split
  minLoss <- function(x) {
    x <- unlist(x)
    if (length(x) > 0) min(x, na.rm = TRUE) else Inf
  }
  vid <- (1:ncol(partvar))[which.min(sapply(loss, minLoss))]
  pid <- (1:nlevels(Part))[which.min(sapply(loss[[vid]], minLoss))]
  pidLab <- nodeids(as.partynode(nodes), terminal = TRUE)[pid]
  stat <- splits[[vid]][[pid]]
  stat <- stat[which.min(stat[, ncol(stat)]), ]
  cut <- stat[1:(length(stat) -1)]
  stat <- stat[length(stat)]
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
  
  ## !!! this part is really ugly and may be improved soon
  
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
               info = list(ordered = ordered,
                 fluctest = test,
                 statistic = stat,
                 splits = splits))
  newnodes[[subs]]$kids <- pidLab + 1L:2L
  newnodes[[subs]]$info$dims <-
    c(N = length(unique(model@subject[subscripts])),
      n = sum(subscripts))
  newnodes[[subs]]$info$step <- step
  
  ## add new children
  newnodes[[subs + 1L]] <- list(id = pidLab + 1L, info = list(dims = c(N = length(unique(model@subject[subsLeft])), n = sum(subsLeft)), depth = newnodes[[subs]]$info$depth + 1L))
  newnodes[[subs + 2L]] <- list(id = pidLab + 2L, info = list(dims = c(N = length(unique(model@subject[subsRight])), n = sum(subsRight)), depth = newnodes[[subs]]$info$depth + 1L))
  
  ## adjust ids
  for (i in 1L:length(newnodes))
    newnodes[[i]]$id <- i
  
  ## print split
  if (control$verbose) {
    
    cat("\n\nSplit = ")
    cat(paste("{",paste(character_split(newnodes[[subs]]$split, data = partvar)$levels, collapse = "}, {"), "}\n", sep = ""))
  }
  
  ## return new nodes
  rval <- as.partynode(newnodes)
  return(rval)
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
  atEtaVar <- terms(model, "fixefEtaVar")
  atEtaInv <- terms(model, "fixefEtaInv")
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
  
  formula <- vcolmm:::olmm_formula(object$info$formula$tree)$fixefEtaVar 
  factors <- attr(terms(formula), "factors")
  if ("Part" %in% all.vars(formula)) {    
    tmp1 <- paste("Eta", 1:model@dims["nEta"], ":", sep  = "")
    tmp2 <- colnames(factors)[factors["Part", ] > 0]
    rval <- c(rval, paste(rep(tmp1, each = length(tmp2)),
                          rep(tmp2, length(tmp1)), sep = ""))
    
  }
  
  formula <- vcolmm:::olmm_formula(object$info$formula$tree)$fixefEtaInv 
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

tvcolmm_fit_checksplit <- function(split, weights, minsplit)
  (sum(split * weights) < minsplit || sum((1 - split) * weights) < minsplit)


tvcolmm_fit_surrogate <- function(split, partvar, weights, model, control) {

  if (control$condsurrogate == "subject") {

    subject <- model@subject
    weights <- unique(cbind(subject, weights))[, -1]
    
    FUN <- function(x) {
      rval <- x[1]
      x <- na.omit(x)
      if (length(x) == 0) rval[1] <- NA else rval <- sample(x, 1)
      return(rval)
    }

    partvar <- aggregate(partvar, by = list(subject), FUN = FUN)
    partvar <- partvar[, -1, drop = FALSE]
  }
  
  x <- partvar[, split$varid]
  xna <- is.na(x)
  
  weights[xna] <- 0L
  
  inp <- rep(TRUE, ncol(partvar))
  inp[split$varid] <- FALSE

  kidsids <- kidids_split(split, data = partvar)
  kidsids[is.na(kidsids)] <- 1L

  if (sum(inp) > 0) {

    rval <-
      partykit:::.csurr(split = kidsids, data = partvar,
                        inp = inp, weights = as.integer(weights),
                        ctrl = list(maxsurrogate = control$maxsurrogate))
  } else {

    rval <- list()
  }
  
  return(rval)
}


tvcolmm_fit_probsplit <- function(split, model, partvar, weights, control) {

  if (control$probsurrogate) {
        
    kidsids <- kidids_split(split, data = partvar)[weights > 0]    

    ## choose randomly a value of an individual
    if (control$condsurrogate == "subject") {
      FUN <- function(x) {
        rval <- x[1]
        x <- na.omit(x)
        if (length(x) == 0) rval[1] <- NA else rval <- sample(x, 1)
        return(rval)
      }
      kidsids <- tapply(kidsids, as.integer(model@subject[weights > 0]), FUN)
    }

    rval <- prop.table(table(kidsids))
    storage.mode(rval) <- "double"
    return(rval)
    
  } else {

    return(NULL)
  }
}


tvcolmm_prune_node <- function(object, alpha = NULL, depth = NULL,
                               minsplit = NULL, nselect = NULL) {

  stopifnot(class(object)[1] %in% c("tvcolmm", "party", "partynode"))
  if ("partynode" %in% class(object)) {
    rval <- object
  } else {
    rval <- object$node
  }

  control <- object$info$control
  if (!is.null(depth) && depth(rval) > 0L)
    rval <- tvcolmm_prune_depth(rval, 0L, depth)
  if (!is.null(minsplit) && depth(rval) > 0L)
    rval <- tvcolmm_prune_minsplit(rval, minsplit)
  if (!is.null(alpha) && depth(rval) > 0L)
    rval <- tvcolmm_prune_alpha(rval, alpha)
  if (!is.null(nselect) && depth(rval) > 0L)
    rval <- tvcolmm_prune_nselect(rval, nselect)

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
      node$info$terminal <- TRUE
      node$info$fluctest <- node$split$info$fluctest
      node$split <- NULL
    }
  }
  return(node)
}


tvcolmm_prune_minsplit <- function(node, minsplit) {
  
  if (!is.terminal(node)) {
    n <- info_node(node)$dims["n"]
    if (n >= 2 * minsplit) {
      kids <- sapply(node$kids, function(kids) kids$id)
      for (i in 1:length(kids)) {
        node$kids[[i]] <- tvcolmm_prune_minsplit(node$kids[[i]], minsplit)
      }
    } else {
      node$kids <- NULL
      node$surrogates <- NULL
      node$info$terminal <- TRUE
      node$info$fluctest <- node$split$info$fluctest
      node$split <- NULL
    }
  }
  return(node)
}


tvcolmm_prune_alpha <- function(node, alpha) {
  
  FUN <- function(node, step) {   
    if (!is.terminal(node)) {
      s <- info_node(node)$step
      if (s < step) {
        kids <- sapply(node$kids, function(kids) kids$id)
        for (i in 1:length(kids)) {
          node$kids[[i]] <- FUN(node$kids[[i]], step)
        }
      } else {
        node$kids <- NULL
        node$surrogates <- NULL
        node$split <- NULL
      }
    }
    return(node)
  }
  
  ids <- setdiff(nodeids(node, terminal = FALSE),
                 nodeids(node, terminal = TRUE))
  steps <- unlist(nodeapply(node, ids, function(node) node$info$step))    
  pval <- extract.tvcolmm(node, "p.value", ids)
  
  if (max(pval) > alpha) {
    step <- min(steps[pval > alpha])
    node <- FUN(node, step)
  } 
  return(node)
}

tvcolmm_prune_nselect <- function(node, nselect) {

  ids <- setdiff(nodeids(node, terminal = FALSE),
                 nodeids(node, terminal = TRUE))    
  pvals <- c(sort(extract.tvcolmm(node, "p.value", ids), decreasing = TRUE), 0)
  selected <- extract.tvcolmm(node, "selected")
  
  while (length(selected) > nselect) {
    node <- tvcolmm_prune_alpha(node, mean(pvals[1:2]))
    ids <- setdiff(nodeids(node, terminal = FALSE),
                   nodeids(node, terminal = TRUE))
    pvals <- c(if (length(ids) > 0) sort(extract.tvcolmm(node, "p.value", ids), decreasing = TRUE), 0)
    selected <- extract.tvcolmm(node, "selected")
  }
  return(node)
}
