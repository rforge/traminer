## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2014-05-02
##
## Description:
## Workhorse functions for the 'tvcm' function
##
## Overview:
##
## tvcm_fit_model:        fits the current model
## tvcm_refit_model:      refits the model
## tvcm_fit_sctest:       run coefficient constancy tests
## tvcm_fit_risk:
## tvcm_fit_splitnode:    split in variable x.
## tvcm_formula:          extract separate formulas for
##                        model and partitioning from
##                        input formula.
## tvcm_update_control:   update the control argument before
##                        fitting the tree.
## tvcm_get_start:        get start values during partitioning
## tvcm_get_node:         extract selected variables
## tvcm_fit_getlevels:
## tvcm_get_terms:       
## tvcm_get_estimates
## tvcm_prune_depth:
## tvcm_prune_minsplit:
## tvcm_prune_minbucket:
## tvcm_prune_maxstep:
## tvcm_prune_terminal:
## tvcm_update_splitpath:
## tvcm_get_estimates:
## tvcm_update_splitpath:
##
## Last modifications:
## 2014-04-27:   complete revision and improved documentation
## 2014-04-01:   rename 'fluctest' to 'sctest'
## 2013-12-02:   remove 'tvcm_fit_setupnode'
## 2013-11-01:   modify 'restricted' and 'terms' correctly in
##               'tvcm_modify_modargs'
##
## To do:
## - improve 'tvcm_refit_model'
## --------------------------------------------------------- #

tvcm_fit_model <- function(call, control) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Fits the current node model.
  ##
  ## Arguments:
  ## call:
  ## control: an object of class 'tvcm_control'.
  ##
  ## Value:
  ## A list of formulas ('root', 'tree' and 'original').
  ##
  ## Details:
  ## Used in 'tvcm' and 'tvcm_fit_sctest'.
  ## ------------------------------------------------------- #
    
  ## extract information from 'call'
  env <- environment(call)
  
  ## fit model
  object <- suppressWarnings(try(eval(call, env), TRUE))

  ## return error if fitting failed
  if (inherits(object, "try-error")) stop("model fitting failed.")

  ## print model
  if (control$verbose) {
    cat("\n\nVarying-coefficient(s) of current model:\n")
    terms <- tvcm_get_terms(names = names(coef(object)),
                            ids = levels(eval(call$data, env)$Node),
                            parm = control$parm)
    print(data.frame(Estimate = coef(object)[terms$node != ""]), digits = 2)
  }
  
  ## return model
  return(object)
}

tvcm_refit_model <- function(object, call) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Re-fits the current node model. Used for the grid-search
  ## in 'tvcm_fit_risk'.
  ##
  ## Arguments:
  ## object:
  ## call:
  ##
  ## Value:
  ## A list of formulas ('root', 'tree' and 'original').
  ##
  ## Details:
  ## Used in 'tvcm' and 'tvcm_fit_risk'. Note that the
  ## function will modify the slots of the original object
  ## as well!
  ##
  ## To do:
  ## Improve performance for non 'olmm' objects
  ## ------------------------------------------------------- #

  env <- environment(call)
  
  if (inherits(object, "olmm")) {
  
    ## set new partition
    slot(object, "frame")$Node <- eval(call$data, env)$Node

    ## get terms
    termsFeCe <- terms(object, "fe-ce")
    termsFeGe <- terms(object, "fe-ge")

    ## get constrasts
    contrasts <- slot(object, "contrasts")
    conCe <- contrasts[intersect(names(contrasts), all.vars(termsFeCe))]
    conGe <- contrasts[intersect(names(contrasts), all.vars(termsFeGe))]    
    
    ## update model matrix  
    slot(object, "X") <-
      olmm_merge_mm(x = model.matrix(termsFeCe, slot(object, "frame"), conCe),
                    y = model.matrix(termsFeGe, slot(object, "frame"), conGe), TRUE)
    
    ## prepare optimization
    optim <- slot(object, "optim")
    optim[[1L]] <- slot(object, "coefficients")
    optim[[4L]] <- slot(object, "restricted")
    environment(optim[[2L]]) <- environment()
    if (!slot(object, "dims")["numGrad"]) environment(optim[[3L]]) <- environment()
    FUN <- optim$fit
    subs <- which(names(optim) == "fit")
    optim <- optim[-subs]
    
    ## run optimization
    output <- suppressWarnings(try(do.call(FUN, optim)))
    if (!inherits(output, "try-error")) {
        slot(object, "output") <- output
    } else {
        object <- output
    }
  } else {
    ## run optimization
    object <- suppressWarnings(try(eval(call, env)))
  }
  
  ## return model
  return(object)
}

tvcm_fit_gefp <- gefp.olmm # for aesthetical reasons ...

tvcm_fit_sctest <- function(model, nodes, partvar, control, call) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Processing of nodewise coefficient constancy tests.
  ##
  ## Arguments:
  ## model:   the current model.
  ## nodes:   an object of class 'partynode'.
  ## partvar: a 'data.frame' with the partitioning variables.
  ## control: an object of class 'tvcm_control'.
  ## call:
  ## Value:
  ## A list with components 'statistic' and 'p.value'.
  ##
  ## Details:
  ## Used in 'tvcm'.
  ## ------------------------------------------------------- #
  
  Node <- factor(fitted_node(nodes, partvar))
  w <- weights(model)
  
  ## get variable types
  functional <- sapply(partvar, function(x) {
    switch(class(x)[1],
           factor = control$functional.factor,
           ordered = control$functional.ordered,
           integer = control$functional.numeric,
           numeric = control$functional.numeric)
  })
  
  
  ## prepare list with arguments for 'sctest'
  dn <- list(paste("Node", levels(Node), sep = ""), colnames(partvar))
  rval <- list(statistic = matrix(, nlevels(Node), ncol(partvar),
                 dimnames = dn),
               p.value = matrix(, nlevels(Node), ncol(partvar),
                 dimnames = dn),
               method = "M-fluctuation test")

  ## get depths of nodes for the check
  depth <- unlist(nodeapply(nodes, nodeids(nodes, terminal = TRUE),
                            function(node) info_node(node)$depth))

  ## extract transformed scores
  eCall <- call(ifelse(inherits(model, "olmm"),"estfun.olmm", "estfun"))
  eCall$x <- quote(model)
  for (arg in names(control$estfun)) eCall[[arg]] <- control$estfun[[arg]]
  scores <- eval(eCall)
  
  ## hack^1: if intercept = "ge" fit a second model with the second level of
  ## 'Node' as reference level
  if (nlevels(Node) > 1L && control$intercept == "ge") {
    contrasts <- eval(call$contrasts, environment(call))
    contrasts$Node <- contr.treatment(levels(Node), base = 2)
    call$contrasts <- contrasts
    verbose <- control$verbose; control$verbose <- FALSE;
    mTmp <- tvcm_fit_model(call, control)
    control$verbose <- verbose
    eCall$x <- quote(mTmp)
    sTmp <- eval(eCall)
    subs1 <- colnames(scores) == paste("Node", levels(Node)[2], sep = "")
    subs2 <- colnames(sTmp) == paste("Node", levels(Node)[1], sep = "")
    attr <- attributes(scores)
    scores <- cbind(scores[, 1:(which(subs1) - 1) ,drop = FALSE],
                    sTmp[, which(subs2) ,drop = FALSE],
                    scores[, which(subs1):ncol(scores) ,drop = FALSE])
    attributes(scores) <- appendDefArgs(attributes(scores), attr)
    rm(mTmp, sTmp)
  }

  terms <- tvcm_get_terms(names = colnames(scores),
                          ids = nodeids(nodes, terminal = TRUE),
                          parm = control$parm)  

  gCall <- call(name = "tvcm_fit_gefp", object = quote(model),
                scores = quote(scores), center = TRUE, silent = TRUE)
  
  ## apply test for each variable and partition separately  
  for (i in 1:ncol(partvar)) {    
    for (j in 1:nlevels(Node)) {
      
      ## extract observations of current partition
      rows <- Node == levels(Node)[j]
        
      ## extract variable to test
      z <- partvar[, i]
      
      ## check if partition contains at least one permitted split
      ## this ensures that we don't select a variable without permitted split
      if (is.numeric(z)) {
        sz <- cumsum(tapply(w[rows], z[rows], sum))
        mb <- max(c(0, apply(cbind(sz, sum(w[rows]) - sz), 1, min)))
      } else {
        sz <- tvcm_fit_getlevels(z, rows)
        mb <- max(c(0, apply(sz, 1, function(x) {
          subs <- z %in% levels(z)[x]
          min(sum(w[subs & rows]), sum(w[!subs & rows]))
        })))
      }
      if (sum(w[rows]) >= control$minsplit & 
          mb >= control$minbucket &
          depth[j] < control$maxdepth) {
      
        ## function to extract scores
        cols <- colnames(scores)[terms$node == levels(Node)[j]]      
  
        ## run the tests
        gCall$order.by <- quote(z)
        gCall$parm <- quote(cols)
        gCall$subset <- quote(rows)
        gefp <- try(eval(gCall), TRUE)
                                    
        if (!inherits(gefp, "try-error")) {

          order.by <- z[rows]

          if (is.character(functional)) {
            functional <- tolower(functional)
            fi <- switch(functional[i],
                         "suplm" = supLM(from = control$trim),
                         "lmuo" = catL2BB(gefp),
                         stop("Unknown efp functional."))
          } else {
            fi <- functional[i]
          }
          test <- try(sctest(x = gefp, functional = fi), TRUE)
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
    for (i in 1:nlevels(Node)) {
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
    if (!is.null(attr(scores, "conv.T")) && attr(scores, "conv.T") == 0L)
      cat("\nNote: Computation of transformation matrix failed.",
          "\nTests are based on raw scores and may be inaccurate.\n")

    if (!is.null(attr(scores, "conv.T")) && attr(scores, "conv.complete") == 0L)
      cat("\nNote: Imputation method for complete data failed.",
          "\nThe p-values may be different in a replication.\n")
    
  }
  return(rval)
}


tvcm_fit_risk <- function(varid, nodeid, model, nodes, partvar, 
                             control, call, step) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Computes the risk for each possible split.
  ##
  ## Arguments:
  ## varid:   an integer vector indicating the partitioning
  ##          variables to be evaluated.
  ## nodeid:  an integer vector indicating the nodes to be
  ##          evaluated.
  ## model:   the current model
  ## nodes:   an object of class 'partynode'.
  ## partvar: a 'data.frame' with the partitioning variables.
  ## control: an object of class 'tvcm_control'.
  ## call:
  ## step:    integer. The current step number.
  ##
  ## Value:
  ## A nested list with risk matrices. Components of nodes
  ## are nested in components for variables. 
  ##
  ## Details:
  ## Used in 'tvcm'.
  ## ------------------------------------------------------- #
  
  if (control$verbose) cat("\nPerform split? ")

  ## if no splitting is necessary ...
  if (length(nodeid) == 0L) {
    if (control$verbose) cat("No\nToo few observations in nodes or/and",
                             "maximum depth reached. Return object.\n")
    return(NULL)
  }

  ## ... otherwise start splitting
  if (control$verbose) cat("Yes\n\n* computing splits ")
  
  Node <- factor(fitted_node(nodes, partvar))    
  w <- weights(model)
  data <- eval(call$data, environment(call))
  call$data <- data
  
  ## function to get all cutpoints of a variable in a partition
  getCutpoints <- function(vid, pid) {
    subs <- Node == levels(Node)[pid]
    if (sum(subs) < 1L | sum(w[subs]) < control$minsplit)
      return(NULL)
    z <- partvar[, vid]
    if (is.numeric(z)) {
        sz <- z[subs]
        w <- w[subs]
        subsL <- which(cumsum(w[order(sz)]) >= control$minbucket)[1]
        subsR <- which(cumsum(w[order(sz)]) > (sum(w) - control$minbucket) - 1)[1]

        if (subsL <= subsR) {            
            rval <- unique(sort(sz)[subsL:subsR])
            if ((length(rval) - 1)  > control$maxevalsplit) 
                rval <-  as.double(quantile(rval, ppoints(control$maxevalsplit)))
            if (rval[length(rval)] == max(z[subs]))
                rval <- rval[-length(rval)]            
        } else {
            rval <- numeric()
        }
        rval <- matrix(rval, ncol = 1)
        colnames(rval) <- "cut";
        if (nrow(rval) > 0L) rownames(rval) <- 1:nrow(rval);
    } else {
        rval <- tvcm_fit_getlevels(z, subs)
    }
    if (nrow(rval) == 0L) rval <- NULL
    return(rval)
  }

  ## function to get the splitstatistic for a given cutpoint
  getSplitstat <- function(cutpoint, vid, pid, call, model) {
    if (!(vid %in% varid & pid %in% nodeid)) return(Inf)
    subs <- Node == levels(Node)[pid]
    z <- partvar[, vid]    
    if (is.numeric(z)) {
      zs <- z <= cutpoint
    } else {
      zs <- z %in% levels(z)[cutpoint]            
    }
    call$data$Node[subs] <- levels(Node)[pid]
    call$data$Node[zs & subs] <- "Right"
    if (sum(w[zs & subs]) < control$minbucket) return(Inf)
    if (sum(w[!zs & subs]) < control$minbucket) return(Inf)
    model <- tvcm_refit_model(model, call)
    rval <- ifelse(inherits(model, "try-error"), NaN, control$riskfun(model))
    if (control$verbose)
      cat(if (is.nan(rval)) "f" else if (rval == Inf) "i" else ".")
    return(rval)
  }

  ## setup a prototyp model to use "tvcm_refit_model" function
  call$data$Node <- Node
  levels(call$data$Node) <- c(levels(call$data$Node), "Right")
  subs <- call$data$Node == levels(call$data$Node)[nodeid[1]]
  call$data$Node[sample(which(subs), as.integer(sum(subs)/2))] <- "Right"
  if (inherits(vcrpart_slot(model, "family"), "family.olmm")) call$doFit <- FALSE
  verbose <- control$verbose; control$verbose <- FALSE;
  splitModel <- tvcm_fit_model(call, control)
  control$verbose <- verbose
  
  ## set start values
  if (inherits(family, "family.olmm") && control$fast > 0L) {
    subs <- intersect(names(coef(splitModel)), names(coef(model)))
    slot(splitModel, "coefficients")[subs] <- slot(model, "coefficients")[subs]
  }
    
  ## run computation (slow!!!)
  riskgrid <- lapply(1:ncol(partvar), function(vid) {
    lapply(1:nlevels(Node), function(pid, vid) {
      cp <- getCutpoints(vid, pid)
      if (is.null(cp)) return(NULL)
      st <- apply(cp, 1, getSplitstat, vid = vid, pid = pid,
                  call, model = splitModel)
      return(cbind(cp, risk = st))
    }, vid = vid)
  })

  Node <- factor(fitted_node(nodes, partvar))
  risk <- lapply(riskgrid, lapply, function(x) x[, ncol(x)])
  
  ## get split
  minRisk <- function(x) {
    x <- unlist(x)
    if (length(x) == 0L) return(Inf)
    x <- na.omit(x)
    if (length(x) > 0L) return(min(x, na.rm = TRUE)) else return(Inf)
  }
  vid <- (1:ncol(partvar))[which.min(sapply(risk, minRisk))]
  pid <- (1:nlevels(Node))[which.min(sapply(risk[[vid]], minRisk))]
  stat <- riskgrid[[vid]][[pid]]
  cid <- which.min(stat[, ncol(stat)])
  
  if (control$verbose) cat(" OK")
  
  return(list(varid = vid, nodeid = pid, cutid = cid, riskgrid = riskgrid))
}

tvcm_fit_splitnode <- function(nodes, risk, partvar, step) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Incorporates a new binary split into an existing
  ## tree structire.
  ##
  ## Arguments:
  ## nodes:   an object of class 'partynode'.
  ## risk:    a list produced by 'tvcm_fit_risk'.
  ## partvar: a 'data.frame' with the partitioning variables.
  ## step:    integer. The current algorithm step.
  ##
  ## Value:
  ## A list of formulas ('root', 'tree' and 'original').
  ##
  ## Details:
  ## Used in 'tvcm'.
  ## ------------------------------------------------------- #
  
  Node <- factor(fitted_node(nodes, partvar))
  vid <- risk$varid
  pid <- risk$nodeid
  pidLab <- nodeids(as.partynode(nodes), terminal = TRUE)[pid]
  stat <- risk$riskgrid[[vid]][[pid]]
  stat <- stat[risk$cutid, ]
  cut <- stat[1:(length(stat) - 1)]
  x <- partvar[, vid]
  
  ## collect information for the split
  subs <- Node == levels(Node)[pid]
  if (is.numeric(x)) { # numerical variables
    breaks <- as.double(cut)
    index <- NULL
    ordered <- TRUE
    subsL <- subs & x <= breaks
    subsR <- subs & x > breaks
  } else {
    subsL <- subs & x %in% levels(x)[cut == 1L]
    subsR <- subs & x %in% levels(x)[cut == 0L]
    if (is.ordered(x)) { 
      breaks <- as.double(max(which(cut == 1)))
      index <- NULL
      ordered <- TRUE
    } else {
      breaks <- NULL
      index <- as.integer(-cut + 2)
      index[table(x[subs]) == 0] <- NA
      ordered <- FALSE
    }
  }
    
  ## get current nodes
  nodes <- as.list(nodes)
  
  ## setup 'newnodes' object
  subsN <- which(sapply(nodes, function(node) node$id) == pidLab)
  newnodes <- vector("list", length(nodes) + 2)
  newnodes[1:subsN] <- nodes[1:subsN]
  if (length(nodes) > pidLab)
    newnodes[(subsN + 3L):length(newnodes)] <-
      nodes[(subsN + 1L):length(nodes)]
  
  ## adjust ids of children
  ids <- sapply(newnodes, function(x) if (!is.null(x$id)) x$id else -1)
  for (i in 1L:length(newnodes))
    if (!is.null(newnodes[[i]]$kids))
      newnodes[[i]]$kids <- which(ids %in% newnodes[[i]]$kids)
  
  ## setup new split
  newnodes[[subsN]]$split <-
    partysplit(varid = vid, breaks = breaks, index = index,
               info = list(ordered = ordered, step = step))
  newnodes[[subsN]]$kids <- pidLab + 1L:2L
  newnodes[[subsN]]$info$dims <-
    c(n = sum(subs))
  
  ## add new children
  newnodes[[subsN + 1L]] <-
    list(id = pidLab + 1L,
         info = list(
           dims = c(n = sum(subsL)),
           depth = newnodes[[subsN]]$info$depth + 1L))
  
  newnodes[[subsN + 2L]] <-
    list(id = pidLab + 2L,
         info =
         list(dims = c(n = sum(subsR)),
              depth = newnodes[[subsN]]$info$depth + 1L))
  
  ## adjust ids
  for (i in 1L:length(newnodes))
    newnodes[[i]]$id <- i
  
  ## return new nodes
  return(as.partynode(newnodes))
}

tvcm_formula <- function(formList, family = cumulative(),
                            env = parent.frame()) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Extracts the formula for the root node and the tree
  ## from the output of 'vcrpart_formula'.
  ##
  ## Arguments:
  ## formList: a list of formulas from 'vcrpart_formula'.
  ## family:   an object of class 'family' or 'family.olmm'.
  ## env:      the environment where the output formula
  ##           is to be evaluated.
  ##
  ## Value:
  ## A list of formulas ('root', 'tree' and 'original').
  ##
  ## Details:
  ## Used in 'tvcm'.
  ## ------------------------------------------------------- #
  
  yName <- all.vars(formList$original)[1] # name of response variable

  ## puts the predictors for fixed effects and varying effects
  ## into one formula
  getTerms <- function(x, type, root = FALSE, family) {
    feTerms <- vcTerms <- NULL
    if (!is.null(x$fe)) {
      fFe <- x$fe$eta[[type]] # formulas for fixed effects
      feTerms <- attr(terms(fFe, keep.order = TRUE), "term.labels")
    }
    if (!is.null(x$vc)) {
      fVc <- x$vc$eta[[type]] # formulas for varying coefficients
      vcTerms <- attr(terms(fVc, keep.order = TRUE), "term.labels")
      if (root) { # delete Node terms
        vcTerms <- vcTerms[vcTerms != "Node"]
        vcTerms <- gsub("Node:", "", vcTerms, fixed = TRUE)
      }
      if ("Node" %in% vcTerms) {
          feTerms <- c("Node", feTerms)
          vcTerms <- vcTerms[vcTerms != "Node"]
      }   
    }
    rval <- ""
    if (length(feTerms) > 0L)
      rval <- paste(rval, paste(feTerms, collapse = "+"), sep = "")
    if (length(feTerms) > 0L & length(vcTerms) > 0L)
      rval <- paste(rval, "+", sep = "")
    if (length(vcTerms) > 0L)
      rval <- paste(rval, paste(vcTerms, collapse = "+"), sep = "")
    if (rval != "" && inherits(family, "family.olmm"))
      rval <- paste(type, "(", rval, ")", sep = "")
    return(rval)
  }

  ## pastes the formulas together
  createFormula <- function(x, root = FALSE, family, env) {
    
    fTree <- "" # the return value

    ## incorporate fixed effects (including varying coefficients)
    feCeTerm <- getTerms(x, "ce", root, family)
    feGeTerm <- getTerms(x, "ge", root, family)
    if (feCeTerm != "") fTree <- paste(fTree, feCeTerm, sep = "")
    if (feCeTerm != "" & feGeTerm != "") fTree <- paste(fTree, " + ", sep = "")
    if (feGeTerm != "") fTree <- paste(fTree, feGeTerm, sep = "")
    if (fTree == "") fTree <- "1"

    ## incorporate intercept
    if (is.null(x$fe) & is.null(x$vc)) {
      x$fe$intercept <- "ce"
    } else if (is.null(x$fe)) {
      x$fe$intercept <- x$vc$intercept
    } else if (is.null(x$vc)) {
      x$vc$intercept <- "none"
    } 
    if (root | x$vc$intercept == "none") {
      if (inherits(family, "family.olmm")) {
        fTree <- paste(fTree, ", intercept='", x$fe$intercept, "'", sep = "")
      } else {
        if (x$fe$intercept == "none")
          fTree <- paste("-1", fTree, sep = "+")
      } 
    } else {
      if (x$vc$intercept == "ce") {
        if (inherits(family, "family.olmm")) {
          fTree <- paste(fTree, ", intercept='none'")
        } else {
          fTree <- paste("-1", fTree, sep = "+")
        }
      }
    }
    if (inherits(family, "family.olmm"))
      fTree <- paste("fe(", fTree, ")", sep = "")

    ## incorporate random effects
    if (!is.null(x$re)) {
      subjectName <- attr(terms(x$re$cond), "term.labels")      
      getReTerms <- function(type) {
        rval <- attr(terms(x$re$eta[[type]], keep.order = TRUE), "term.labels")
        if (x$re$intercept == type) rval <- c("1", rval)
        if (length(rval) == 0L) return(NULL)
        rval <- paste(rval, collapse = "+")
        if (inherits(family, "family.olmm"))
          rval <- paste(type, "(", rval, ")", sep = "")
        return(rval)
      }
      reTerm <- unlist(lapply(c("ce", "ge"), getReTerms))
      reTerm <- paste(reTerm, collapse = "+")
      if (inherits(family, "family.olmm")) {        
        reTerm <- paste("re(", reTerm, "|", subjectName,
                        ",intercept='", x$re$intercept, "')", sep = "")
      } else {        
        if (x$re$intercept == "none")
          reTerm <- paste(reTerm, "-1", sep = "")
        reTerm <- paste("(", reTerm, "|", subjectName, ")", sep = "")
      }        
      fTree <- paste(fTree, "+", reTerm)
    }

    ## incorporate response variable and create a formula
    fTree <- paste(yName, "~", fTree)
    return(formula(fTree))
  }

  rval <- list(root = createFormula(formList, root = TRUE, family, env),
               tree = createFormula(formList, root = FALSE, family, env),
               original = formList$original)
  environment(rval$root) <- environment(rval$tree) <- env
  return(rval)
}

tvcm_update_control <- function(control, model, formList) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Adds a new slot 'parm' to the 'control_tvcm' object
  ## for internal purposes.
  ##
  ## Arguments:
  ## control:  an object of class 'tvcm_control'.
  ## model:    a root node regression model, e.g., an 'olmm'
  ##           or a 'glm' object
  ## formList: a list of formulas from 'vcrpart_formula'.
  ## 
  ##
  ## Value:
  ## An updated 'tvcm_control' object.
  ##
  ## Details:
  ## Used in 'tvcm'.
  ## ------------------------------------------------------- #

  family <- vcrpart_slot(model, "family")
  if (is.null(formList$vc)) return(control)

  ## set 'vcterms' slot
  vcTerms <- list(ce = attr(terms(formList$vc$eta$ce), "term.labels"),
                  ge = attr(terms(formList$vc$eta$ge), "term.labels"))
  getNewTerms <- function(term, type) {
    if (term == "Node") return(NULL)
    terms <- terms(model, which = type)
    X <- model.matrix(model, which = type)
    assign <- attr(X, "assign")
    term <- gsub("Node:", "", term)
    subs <- which(attr(terms, "term.label") == term)
    if (length(subs) < 1L) return(NULL)
    rval <- colnames(X)[assign == subs]
    return(rval)
  }
  vcTerms <- lapply(names(vcTerms), function(i) if (length(vcTerms[[i]]) > 0L) return(unlist(lapply(1:length(vcTerms[[i]]), function(j) getNewTerms(vcTerms[[i]][j], paste("fe", i, sep = "-"))))) else return(NULL))
  if (inherits(family, "family.olmm") & length(vcTerms[[1L]]) > 0L)
    vcTerms[[1L]] <- paste("Eta", 1L:slot(model, "dims")["nEta"], ":", vcTerms[[1L]], sep = "")
  control$parm <- unlist(vcTerms)
  
  ## set 'intercept' slot
  control$intercept <- formList$vc$intercept
  if (control$intercept != "none") {
    if (inherits(model, "olmm")) {
      if (control$intercept == "ce") {
        control$parm <- c(grep("Eta[1-9]+:\\(Intercept\\)",
                               names(coef(model)), value = TRUE),
                          control$parm)
      }
    } else {
      control$parm <- c("(Intercept)", control$parm)
    }
  }

  ## set 'nuisance' slot for estfun
  control$estfun$nuisance <- setdiff(names(coef(model)), control$parm)
  return(control)
}

tvcm_get_start <- function(model, ids, parm) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Extracts non-varying coefficients of a current 'model'
  ## to be used as intial values for re-fitting the model.
  ##
  ## Arguments:
  ## model: the current model
  ## ids:
  ## parm:
  ## 
  ## Value:
  ## A named numeric vector with initial values for fixed
  ## effects excluding varying effects.
  ##
  ## Details:
  ## Used in 'tvcm'.
  ## ------------------------------------------------------- #
  
  coef <- coef(model)
  terms <- tvcm_get_terms(names(coef), ids, parm)
  return(coef[terms$type != "vc"])
}

tvcm_get_node <- function(object, newdata, setContrasts = FALSE, weights) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Extract the node vector from 'newdata' and assign
  ## the contrasts.
  ##
  ## Arguments:
  ## object:  a fitted 'tvcm' object.
  ## newdata: a data.frame from which the nodes are to be
  ##          extracted.
  ## weights: a numeric vector of weights with the same
  ##          size than data.
  ##
  ## Value:
  ## A list with values of the variables in 'data'
  ##
  ## Details:
  ## Used in tvcm.predict and tvcm.prune.
  ## ------------------------------------------------------- #  
  
  fitted <-
    fitted_node(object$node, newdata[, colnames(object$data), drop = FALSE])
  rval <- factor(fitted, levels = nodeids(object$node, terminal = TRUE))
  names(rval) <- rownames(newdata)

  if (nlevels(rval) > 1L) {
    if (setContrasts) {
      if (is.null(weights)) stop("weights are required")
      con <- contr.sum(levels(rval))
      tab <- tapply(weights, rval, sum)
      con[nrow(con),] <- con[nrow(con),] * tab[-length(tab)] / tab[length(tab)]
      colnames(con) <- levels(rval)[1:(nlevels(rval) - 1)]
      contrasts(rval) <- con
    } else {
      contrasts(rval) <- vcrpart_slot(object, "contrasts")$Node
    }
  }
  return(rval)
}

tvcm_fit_getlevels <- function(x, subset) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Computes all possible binary splits for a categorical
  ## variable
  ##
  ## Arguments:
  ## x:  
  ## subset: a vector of subscripts for the subset for
  ##             which the splits should be computed
  ##
  ## Value:
  ## A matrix with splits
  ##
  ## Details:
  ## Used in tvcm_fit_sctest and tvcm_fit_risk.
  ## Function is a modified version of the 'mob_fit_getlevels'
  ## function of the package 'party'.
  ## ------------------------------------------------------- #
  
  ## all levels
  nl <- nlevels(x)

  ## levels in current partition
  xd <- droplevels(x[subset > 0])
  nld <- nlevels(xd)
  
  if (inherits(x, "ordered")) {
    indx <- diag(nl)
    indx[lower.tri(indx)] <- 1
    indx <- indx[-nl, , drop = FALSE]
    rownames(indx) <- levels(x)[-nl]
  } else {
    mi <- 2^(nld - 1) - 1
    indx <- matrix(0, nrow = mi, ncol = nl)
    for (i in seq(1, mi, length.out = mi)) { # go through all splits
      ii <- i
      for (l in 1:nld) {
        indx[i, which(levels(x) %in% levels(xd))[l]] <- ii %% 2;
        ii <- ii %/% 2   
      }
    }
    rownames(indx) <-
      apply(indx, 1, function(z) paste(levels(x)[!is.na(z) & z > 0], collapse = "+"))
  }
  colnames(indx) <- as.character(levels(x))
  storage.mode(indx) <- "logical"
  return(indx)
}

tvcm_get_terms <- function(names, ids, parm) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Creates a list which can be used to extract the varying
  ## coefficients.
  ##
  ## Arguments:
  ## names: character vector. Names of coefficients of the
  ##        current model.
  ## ids:   character vector. Names of the current nodes.
  ##
  ## Value:
  ## A list with slots
  ## names:  the original coefficient names.
  ## terms:  the names of the terms to which coefficients
  ##         belong according to the original formula.
  ## type:   which type of 'fe', 'vc' or 're' the term
  ##         belongs to.
  ## node:   the node to which a coefficient belongs to.
  ## ------------------------------------------------------- #

  if (length(ids) == 1L && ids == 1L)
      names[names %in% parm] <- paste("Node1", names[names %in% parm], sep = ":")
  nodes <- paste("Node", ids, sep = "")
  split <- strsplit(names, ":")
  type <-
    sapply(split, function(x) {
      rval <- "fe"
      if (any(x %in% nodes)) rval <- "vc"
      if (any(substr(x, 1, 12) %in% "ranefCholFac")) rval <- "re"
      return(rval)
    })
  terms <-
    sapply(split, function(x) {
      paste(x[!x %in% nodes], collapse = ":")
    })
  node <- 
    sapply(split, function(x) {
      rval <- ""
      if (any(subs <- x %in% nodes))
        rval <- substr(x[subs], 5, 500)
      return(rval)
    })
  return(list(names = names, terms = terms, type = type, node = node))
}

tvcm_get_estimates <- function(object, what = c("coef", "sd", "var"), ...) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Extracts the estimates a fitted 'tvcm' object and
  ## creates a list used in various methods.
  ## variances of a fitted 'tvcm' object.
  ##
  ## Arguments:
  ## object: an object of class 'partynode'.
  ## what:   the type of estimate to be extracted.
  ##
  ## Value:
  ## A list of class 'splitpath.tvcm'.
  ##
  ## Details:
  ## Used in 'coef', 'extract'.
  ## ------------------------------------------------------- #
  
  what <- match.arg(what)
  model <- extract(object, "model")

  rval <- list(fe = numeric(), vc = matrix(,0,0), re = numeric())
  
  ## extract coefficients
  estimates <- switch(what,
                      coef = coef(model),
                      sd = diag(vcov(model)),
                      var = diag(vcov(model)))

  ids <- nodeids(object, terminal = TRUE)
  names <- names(estimates)
  terms <- tvcm_get_terms(names, ids, object$info$control$parm)
  
  ## restricted coefficients
  if (any(terms$type == "fe"))
    rval$fe <- estimates[terms$type == "fe"]
  if (any(terms$type == "re"))
    rval$re <- estimates[terms$type == "re"]
  
  ## varying coefficients
  if (any(terms$type == "vc")) {
    vcTerms <- unique(terms$terms[terms$type == "vc"])
    rval$vc <- matrix(, length(ids), length(vcTerms))
    rownames(rval$vc) <- ids; colnames(rval$vc) <- vcTerms;
    for (i in 1:ncol(rval$vc)) {
      subs <- terms$terms == vcTerms[i] & terms$node != ""
      rval$vc[terms$node[subs], i] <- estimates[subs]
    }
    subs <- colnames(rval$vc) %in% c("")
    colnames(rval$vc)[subs] <- "(Intercept)"
    if (inherits(object$info$family, "family.olmm")) {
      subs <- substr(colnames(rval$vc), 1L, 3L) == "Eta"
      colnames(rval$vc)[subs] <-
        paste(colnames(rval$vc)[subs], "(Intercept)", sep = ":")
    }
    if ("(Intercept)" %in% colnames(rval$vc)) {
      con <- vcrpart_slot(model, "contrasts")$Node
      if (is.character(con))
        con <- do.call(con, list(n = ids))
      intercepts <- rval$vc[, "(Intercept)"]
      subs <- is.na(intercepts)
      if (what == "coef") {
        intercepts[subs] <- sum(con[subs, ] * intercepts[!subs])
        rval$vc[, "(Intercept)"] <- intercepts
      } else if (what %in% c("sd", "var")) {
        intercepts[subs] <- sum(con[subs, ]^2 * intercepts[!subs])
        rval$vc[, "(Intercept)"] <- intercepts
      }
    }
    if (what == "sd")
      rval <- lapply(rval, function(x) if (!is.null(x)) sqrt(x) else x)

  } else {
    rval$vc <- matrix(, 0, 0)
  }  
  return(rval)
}

tvcm_prune_node <- function(object, keepids = FALSE, alpha = NULL,
                               maxdepth = NULL, maxwidth = NULL,
                               minsplit = NULL, minbucket = NULL,
                               nselect = NULL, maxstep = NULL,
                               terminal = NULL) {

  ## ------------------------------------------------------- #
  ## Description:
  ## The main function to prune an object of class 'partynode'
  ## according to some criteria.
  ##
  ## Arguments:
  ## nodes:      an object of class 'partynode'.
  ## keepids:    whether the original ids of the unpruned
  ##             'nodes' object should be inscribed in the
  ##             info slots of the according nodes.
  ## alpha, ...: criteria according to which 'nodes' is to be
  ##             pruned.
  ##
  ## Value:
  ## A list of class 'splitpath.tvcm'.
  ##
  ## Details:
  ## Used in 'tvcm'.
  ## ------------------------------------------------------- #
  
  stopifnot(class(object)[1] %in% c("tvcm", "party", "partynode"))
  if ("partynode" %in% class(object)) {
    rval <- object
  } else {
    rval <- object$node
  }

  control <- extract(object, "control")
  if (!is.null(maxdepth) && depth(rval) > 0L && maxdepth >= 0L)
    rval <- tvcm_prune_depth(rval, 0L, maxdepth)
  if (!is.null(maxwidth) && width(rval) > 0L && maxwidth > 0L)
      rval <- tvcm_prune_maxstep(rval, maxwidth - 1L)
  if (!is.null(minsplit) && depth(rval) > 0L)
    rval <- tvcm_prune_minsplit(rval, minsplit)
  if (!is.null(minbucket) && depth(rval) > 0L)
    rval <- tvcm_prune_minbucket(rval, minbucket)
  if (!is.null(nselect) && length(extract(object, "selected")) > 0L) {
    splitpath <- object$info$splitpath
    vars <- sapply(splitpath, function(x) x$varid)
    nvars <- sapply(1:length(vars), function(x) length(unique(vars[1:x])))
    maxstep <- max(c(0, which(nvars <= nselect)))
    rval <- tvcm_prune_maxstep(rval, maxstep)
  }
  if (!is.null(alpha) && depth(rval) > 0L) {
    splitpath <- object$info$splitpath
    p.value <- extract(object, "p.value")
    maxstep <- c(0, which(!is.na(p.value) & p.value <= alpha))
    if (length(maxstep) > 0L) 
      maxstep <- max(maxstep[c(1, diff(maxstep)) == 1])
    rval <- tvcm_prune_maxstep(rval, maxstep)
  }
  if (!is.null(maxstep) && depth(rval) > 0L) {
    rval <- tvcm_prune_maxstep(rval, maxstep)
  }
  if (!is.null(terminal) && depth(rval) > 0L) {
    rval <- tvcm_prune_terminal(rval, terminal)
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
  if (!keepids) {
      for (i in 1:length(rval)) rval[[i]]$info$original$id <- rval[[i]]$id
      rval <- as.partynode(rval)
  }
  rval <- as.partynode(rval)
  return(rval)
}

tvcm_prune_depth <- function(node, d, maxdepth) {
  
  if (!is.terminal(node)) {
    if (d + 1 <= maxdepth) {
      kids <- sapply(node$kids, function(kids) kids$id)
      for (i in 1:length(kids)) {
        node$kids[[i]] <-
          tvcm_prune_depth(node$kids[[i]], d + 1L, maxdepth)
      }
    } else {
      node$kids <- NULL
      node$split <- NULL
    }
  }
  return(node)
}

tvcm_prune_minsplit <- function(node, minsplit) {
  
  if (!is.terminal(node)) {
    n <- info_node(node)$dims["n"]
    if (n >= minsplit) {
      kids <- sapply(node$kids, function(kids) kids$id)
      for (i in 1:length(kids)) {
        node$kids[[i]] <- tvcm_prune_minsplit(node$kids[[i]], minsplit)
      }
    } else {
      node$kids <- NULL
      node$split <- NULL
    }
  }
  return(node)
}

tvcm_prune_minbucket <- function(node, minbucket) {
  
  if (!is.terminal(node)) {
    n <- sapply(node$kids, function(x) x$info$dims["n"])
    if (min(n) >= minbucket) {
      kids <- sapply(node$kids, function(kids) kids$id)
      for (i in 1:length(kids)) {
        node$kids[[i]] <- tvcm_prune_minbucket(node$kids[[i]], minbucket)
      }
    } else {
      node$kids <- NULL
      node$split <- NULL
    }
  }
  return(node)
}

tvcm_prune_maxstep <- function(node, maxstep) {
  
  if (!is.terminal(node)) {
    if (node$split$info$step <= maxstep) {
      kids <- sapply(node$kids, function(kids) kids$id)
      for (i in 1:length(kids)) {
        node$kids[[i]] <- tvcm_prune_maxstep(node$kids[[i]], maxstep)
      }
    } else {
      node$kids <- NULL
      node$split <- NULL
    }
  }
  return(node)
}

tvcm_prune_terminal <- function(node, terminal) {
  
  if (!is.terminal(node)) {
    if (!node$id %in% terminal) {
      kids <- sapply(node$kids, function(kids) kids$id)
      for (i in 1:length(kids)) {
        node$kids[[i]] <- tvcm_prune_terminal(node$kids[[i]], terminal)
      }
    } else {
      node$kids <- NULL
      node$split <- NULL
    }
  }
  return(node)
}

tvcm_update_splitpath <- function(splitpath, nodes, partvar, control) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Creates a 'splitpath.tvcm' object for tracing the
  ## fitting process.
  ##
  ## Arguments:
  ## splitpath: the splitpath obtained by the fitting process.
  ## nodes:     an object of class 'partynode'.
  ## partvar:   a 'data.frame' with the partitioning variables.
  ## control:   an object of class 'tvcm_control'.
  ##
  ## Value:
  ## A list of class 'splitpath.tvcm'.
  ##
  ## Details:
  ## Used in 'tvcm'.
  ## ------------------------------------------------------- #
  
  if (width(nodes) < 2L) return(splitpath)
  ids <- nodeids(nodes)
  steps <- nodeapply(nodes,  ids, function(node) node$split$info$step)
  steps <- sapply(steps, function(x) if (is.null(x)) Inf else x)
  
  for (step in 1:sum(steps < Inf)) {
    if (step > 1L) {
      parentids <- ids[steps < step]
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
    if (control$method == "mob") {
      if (!is.null(splitpath[[step]]$sctest)) {
        rownames(splitpath[[step]]$sctest$statistic) <-
          paste("Node", kidids, sep = "")
        rownames(splitpath[[step]]$sctest$p.value) <-
          paste("Node", kidids, sep = "")
      }
      if (!is.null(splitpath[[step]]$varid)) {
        if (length(splitpath[[step]]$riskgrid) > 1L)
          splitpath[[step]]$riskgrid[-splitpath[[step]]$varid] <- NULL
        if (length(splitpath[[step]]$riskgrid[[1]]) > 1L)
          splitpath[[step]]$riskgrid[[1]][-splitpath[[step]]$nodeid] <- NULL
        names(splitpath[[step]]$riskgrid)[1L] <- 
          colnames(partvar)[splitpath[[step]]$varid]
        names(splitpath[[step]]$riskgrid[[1]])[1L] <- 
          paste("Node", ids[steps == step], sep = "")
      }
    } else {
      if (!is.null(splitpath[[step]]$riskgrid)) {
        names(splitpath[[step]]$riskgrid) <- colnames(partvar)
        for (i in 1:ncol(partvar)) 
          names(splitpath[[step]]$riskgrid[[i]]) <- 
            paste("Node", kidids, sep = "")
      }
    }
    if (!is.null(splitpath[[step]]$varid)) {
      splitpath[[step]]$nodeid <- ids[steps == step]
      splitpath[[step]]$var <- colnames(partvar)[splitpath[[step]]$varid]
      splitpath[[step]]$cutpoint <- character_split(nodeapply(nodes, ids[steps == step], function(node) node$split)[[1]], data = partvar)$levels
    }
  }
  class(splitpath) <- "splitpath.tvcm"
  return(splitpath)
}
