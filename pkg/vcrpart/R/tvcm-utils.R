##' -------------------------------------------------------- #
##' Author:          Reto Buergin
##' E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
##' Date:            2014-09-08
##'
##' Description:
##' Workhorse functions for the 'tvcm' function
##'
##' Overview:
##'
##' Workhorse functions for partitioning:
##' tvcm_grow:                main function for growing the trees
##' tvcm_grow_fit:            fits the current model
##' tvcm_grow_update:          refits the model (with utility function
##'                           'glm.doNotFit')
##' tvcm_grow_setsplits:      get current splits
##' tvcm_setsplits_validcats: validate categorical cuts
##' tvcm_setsplits_sctest:    update splits after tests
##' tvcm_setsplits_splitnode: 
##' tvcm_setsplits_rselect:   randomly select partitions, variables and nodes
##' tvcm_grow_sctest:         run coefficient constancy tests
##' tvcm_grow_loss:           grid based loss minimization
##' tvcm_grow_splitnode:      split in variable x.
##' tvcm_formula:             extract separate formulas for
##'                           model and partitioning from
##'                           input formula.
##' tvcm_grow_setcontrol:     update the control argument before
##'                           fitting the tree.
##'
##' Utility functions used by various functions:
##' tvcm_get_node:            extract node vectors and assign the contrasts
##' tvcm_get_terms:           creates a list which assigns coefficients
##'                           to the corresponding type, partition etc.
##' tvcm_get_vcparm:          extracts the names of the predictors on 'vc' terms
##' tvcm_get_estimates:       extracts the estimates from a fitted
##'                           'tvcm' object and creates a list with
##'                           an entry for each different type of
##'                           estimate ('fe', 'vc' or 're')
##' tvcm_print_vclabs:        creates short labels for 'vc' terms
##'
##' Functions for pruning:
##' tvcm_prune_node:          main function for pruning 'partynode'
##'                           objects
##' tvcm_prune_maxstep:       recursive function for pruning
##' tvcm_prune_terminal:      prunes branches
##' tvcm_grow_splitpath:      creates a 'splitpath.tvcm' object
##'
##' Last modifications:
##' 2014-09-08: substitute 'rep' function by 'rep.int' or 'rep_len'
##' 2014-09-07: - added 'tvcm_get_vcparm' function
##'             - set default values in 'glm.doNotFit'
##' 2014-09-06: modified function names for 'tvcm_fit_model' and
##'             'tvcm_refit_model' for consistency reasons. The
##'             new names are 'tvcm_grow_fit' and 'tvcm_grow_update'
##' 2014-09-06: added new function 'tvcm_grow', which was formerly
##'             in 'tvcm'
##' 2014-09-04: added new function 'tvcm_print_vclabs'
##' 2014-09-02: modifications on 'tvcm_get_node' to accelerate
##'             the code
##' 2014-08-10: modifications to speed-up the code
##'             - update formulas of 'tvcm_formula' are now
##'               always identical
##'             - if 'doFit = FALSE', the call of 'glm.fit'
##'               is avoided
##' 2014-08-08: correct bug in 'tvcm_get_terms' for cases where
##'             multiple vc() terms with equal 'by' arguments
##'             are present
##' 2014-08-08: correct bug in 'tvcm_grow_setsplits' regarding
##'             'keeploss'
##' 2014-08-08: add suppressWarnings in tvcm_grow_fit
##' 2014-07-22: the list of splits is now of structure
##'             partitions-nodes-variables
##' 2014-07-22: AIC and BIC are no longer criteria and therefore
##'             multiple functions were adjusted
##' 2014-07-22: modified some function names
##' 2014-07-06: implement method to deal with many nominal categories
##' 2014-06-30: implement random selection if split is not unique
##' 2014-06-23: correct bug for 'start' argument in 'tvcm_grow_loss' 
##' 2014-06-17: modify documentation style
##' 2014-06-16: deleted several 'tvcm_prune_XXX' functions
##' 2014-06-03: modify 'tvcm_formula' to allow partition-wise
##'             trees
##' 2014-04-27: complete revision and improved documentation
##' 2014-04-01: rename 'fluctest' to 'sctest'
##' 2013-12-02: remove 'tvcm_grow_setupnode'
##' 2013-11-01: modify 'restricted' and 'terms' correctly in
##'             'tvcm_modify_modargs'
##' -------------------------------------------------------- #

##' -------------------------------------------------------- #
##' \code{\link{tvcm_grow_fit}} fits the current node model.
##'
##' @param object  a 'tvcm' object
##' @param subset  a vector indicating the subset on which
##'    the model is to be fitted
##' @param weights a vector of weights corresponding to the
##'    subset entries
##' 
##' @return A 'tvcm' object.
##'
##' @details Used in 'cvloss.tvcm' 
##' -------------------------------------------------------- #

tvcm_grow <- function(object, subset = NULL, weights = NULL) {

  mcall <- object$info$mcall
  environment(mcall) <- environment()
  formList <- object$info$formula
  model <- object$info$model
  mf <- model.frame(model)
  partData <- object$data
  control <- object$info$control
  family <- model$family
  start <- object$dotargs$start
  for (arg in names(object$info$dotargs))
    assign(arg, object$info$dotargs[[arg]])
  
  if (!is.null(subset)) {
    mf <- mf[subset,,drop = FALSE]
    partData <- partData[subset,, drop = FALSE]
  }

  if (!is.null(weights)) {
    mcall$weights <- weights
  } else {
    weights <- weights(model)
  }
      
  ## get number of partitions
  nPart <- length(formList$vc)

  ## get partitioning variables
  partVars <- lapply(formList$vc, function(x) attr(terms(x$cond), "term.labels"))
  varid <- lapply(partVars, function(x) {
    as.integer(sapply(x, function(x) which(colnames(partData) == x))) })
     
  ## set the root node
  nodes <- replicate(nPart, partynode(id = 1L, info = list(dims = nobs(model), depth = 0L)))
  names(nodes) <- names(formList$vc)
  where <- vector("list", length = nPart)
  
  partid <- seq(1, nPart, length.out = nPart)
  spart <- 0 # pseudo value
  splits <- vector("list", length = nPart)
  
  splitpath <- list()
  
  run <- 1L
  step <- 0L
  noverstep <- 0L
  
  while (run > 0L) {

    step <- step + 1L; nstep <- step;
    test <- NULL; loss <- NULL;

    ## get current partitions and add them to the model data
    for (pid in seq_along(nodes)) {
      where[[pid]] <- factor(fitted_node(nodes[[pid]], partData))
      where[[pid]] <- vcrpart_contr.sum(where[[pid]], weights)      
      mf[, paste("Node", LETTERS[pid], sep = "")] <- where[[pid]]
    }
    
    nodeid <- lapply(nodes, function(x) 1:width(x))
    
    if (control$verbose) cat("\n* starting step", step, "...")
    
    ## --------------------------------------------------- #
    ## Step 1: fit the current model
    ## --------------------------------------------------- #

    vcRoot <- sapply(nodeid, length) == 1L
    ff <- tvcm_formula(formList, vcRoot, model$family,
                       environment(formList$original))
    model <- try(tvcm_grow_fit(mcall))
    
    if (inherits(model, "try-error")) stop(model)

    control <- tvcm_grow_setcontrol(control, model, formList, vcRoot)

    if (control$verbose) {
      cat("\n\nVarying-coefficient(s) of current model:\n")
      if (length(unlist(control$parm)) > 0L) {
        print(data.frame(Estimate = coef(model)[unique(unlist(control$parm))]),
              digits = 2)
      } else {
        cat("<no varying-coefficients>\n")
      }
    }
    
    ## compute / update splits
    splits <- tvcm_grow_setsplits(splits, spart, partid, nodeid, varid, model,
                                  nodes, where, partData, control)
    
    ## check if there is at least one admissible split
    if (length(unlist(splits)) == 0L | step > control$maxstep |
        length(control$parm) == 0L) {
      run <- 0L
      if (step > control$maxstep) {
        stopinfo <- "maximal number of steps reached"
      } else if (length(control$parm) == 0L) {
        stopinfo <- "no varying coefficients"
      } else {
        stopinfo <- "no admissible splits (exceeded tree size parameters)"
      }
      nstep <- step - 1L
    }

    ## random selection (used by 'fvcm')
    if (any(c(control$ptry, control$vtry, control$ntry) < Inf))
        splits <- tvcm_setsplits_rselect(splits, partid, nodeid, varid, control)
    
    if (run > 0L && control$sctest) {
      
      ## --------------------------------------------------- #
      ## Step 2: variable selection via coefficient constancy tests
      ## --------------------------------------------------- #
      
      ## get raw p-values
      test <- try(tvcm_grow_sctest(model, nodes, where, partid, nodeid, varid, 
                                   splits, partData, control), TRUE)
      
      ## return error if test failed
      if (inherits(test, "try-error")) {
        run <- 0L
        stopinfo <- test

      } else {      
        testAdj <- tvcm_sctest_bonf(test,ifelse(control$bonferroni,"nodewise", "none"))
        run <- 1L * (min(c(1.0 + .Machine$double.eps, unlist(testAdj)),
                         na.rm = TRUE) <= control$alpha)
      }
      
      if (run > 0L) {
        
        ## extract the selected partition
        testAdjPart <-
          tvcm_sctest_bonf(test,ifelse(control$bonferroni,"partitionwise","none"))
        minpval <- min(unlist(testAdjPart), na.rm = TRUE)
        spart <- which(sapply(testAdjPart, function(x)any(sapply(x,identical,minpval))))
        if (length(spart) > 1L) spart <- sample(spart, 1L)

        ## select variable and node
        minsubs <- which(sapply(test[[spart]], identical,
                                min(test[[spart]], na.rm = TRUE)))
        if (length(minsubs) > 1L) minsubs <- sample(minsubs, 1L)
        svar <- ceiling(minsubs / nrow(test[[spart]]))
        snode <- minsubs - (svar - 1L) * nrow(test[[spart]])
        
        ## print results
        if (control$verbose) {

          ## tests
          cat("\nCoefficient constancy tests (p-value):\n")   
          for (pid in seq_along(nodes)) {
                cat(paste("\nPartition ", LETTERS[pid], ":\n", sep = ""))
                print(data.frame(format(testAdj[[pid]], digits = 2L)))              
              }
          
          ## selections
          cat("\nSplitting partition:", names(nodes)[spart])
          cat("\nSplitting variable:", names(partData)[varid[[spart]][svar]])
          cat("\nNode:", levels(where[[spart]])[snode])
          
        }

        ## set deviance statistic of not to selected nodes to 'Inf' to avoid
        ## model evaluations
        splits <- tvcm_setsplits_sctest(splits, partid, spart,
                                        nodeid, snode, varid, svar)
        
      } else {
        stopinfo <- "p-values of coefficient constancy tests exceed alpha"
      }  
    }
      
    if (run > 0L) {
      
      ## ------------------------------------------------- #
      ## Step 3: search a cutpoint
      ## ------------------------------------------------- #
      
      ## compute the loss of all candidate splits and extract the best split
      loss <- try(tvcm_grow_loss(splits, partid, nodeid, varid, 
                                 model, nodes, where, partData,
                                 control, mcall, formList, step), silent = TRUE)
      
      ## handling stops
      if (inherits(loss, "try-error")) {
        run <- 0L
        stopinfo <- loss
        nstep <- step - 1L
        
      } else {
        splits <- loss$lossgrid
        spart <- loss$partid
        
        if (is.null(loss$cut)) {
          run <- 0L
          stopinfo <- "no split that decreases the loss found"
          nstep <- step - 1L
        }
        
        if (run > 0L) {
          noverstep <- if (loss$ploss < 0) noverstep + 1L else 0L
          if (noverstep > control$maxoverstep) {
            run <- 0L
            if (control$maxoverstep > 0L) {
                stopinfo <- "'maxoverstep' reached"
            } else {
                stopinfo <- "minimal penalized loss-reduction reached"
            }
            nstep <- nstep - 1L
          }
        }
        
      }
    }
    
    ## incorporate the split into 'nodes'
    if (run > 0L)
      nodes <- tvcm_grow_splitnode(nodes, where, loss, partData,
                                   step, weights)

    if (run > 0L)
      splits <- tvcm_setsplits_splitnode(splits, loss$partid, loss$nodeid,
                                         nodeid, loss, model, control)
      
    ## update 'splitpath' to make the splitting process traceable
    if (run >= 0L)
      splitpath[[step]] <-
        list(step = step,
             loss = control$lossfun(model),
             npar = extractAIC(model)[1L],
             nspl = step - 1L)

    if (!inherits(test, "try-error"))
      splitpath[[step]]$sctest <- test
    
    if (!inherits(loss, "try-error")) {
      if (run > 0L) {
        splitpath[[step]]$partid <- loss$partid
        splitpath[[step]]$nodeid <- loss$nodeid
        splitpath[[step]]$varid <- loss$varid
        splitpath[[step]]$cutid <- loss$cutid
        splitpath[[step]]$ploss <- loss$ploss
      }
      splitpath[[step]]$lossgrid <- loss$lossgrid 
    }
    
    ## print the split
    if (control$verbose) {
      if (run > 0L) {
        if (!control$sctest) {
            cat("\n\nSplitting partition:", names(nodes)[loss$partid])
            cat("\nNode:", levels(where[[loss$partid]])[loss$nodeid])
            cat("\nVariable:", names(partData)[loss$varid])
        } else {
            cat("\n")
        }

        cat("\nCutpoint:\n")
        print(as.data.frame(matrix(loss$cut, 1L,
                                   dimnames = list(loss$cutid,
                                     names(loss$cut)))))
        
        cat("Model comparison:\n")
        print(data.frame("loss" = c(control$lossfun(model),
                           control$lossfun(model) - loss$loss),
                         "penalized loss reduction" =  c("", format(loss$ploss)),
                         row.names = paste("step", step + c(-1, 0)),
                         check.names = FALSE))
        
      } else {
        cat("\n\nStopping the algorithm.\nMessage:", as.character(stopinfo), "\n")
        if (inherits("try-error", stopinfo)) warning(as.character(stopinfo))
        
      }
    }
  }
  
  ## inscribe original node names for later pruning
  for (pid in seq_along(nodes)) {
    nodes[[pid]] <- as.list(nodes[[pid]])
    for (nid in 1:length(nodes[[pid]])) {
      nodes[[pid]][[nid]]$info$id$original <- nodes[[pid]][[nid]]$id
      nodes[[pid]][[nid]]$info$id$last <- nodes[[pid]][[nid]]$id
    }
    nodes[[pid]] <- as.partynode(nodes[[pid]])
  }
  
  ## prepare the title
  title <- c("Tree-based varying-coefficients model")
  
  ## modify splitpath    
  splitpath <- tvcm_grow_splitpath(splitpath, varid, nodes, partData, control)
  
  ## the output object
  if (nPart == 0L) {
    tree <- model
    
  } else {

    ## delete environments of calls
    environment(mcall) <- NULL
    environment(object$info$call) <- NULL

    ## build 'tvcm' object
    tree <- party(nodes[[1L]],
                  data = partData,
                  fitted = data.frame(
                    "(response)" = model.response(model.frame(model)),
                    "(weights)" = weights(model),
                    check.names = FALSE),
                  terms = terms(formList$original, keep.order = TRUE),
                  info = list(
                    title = title,
                    call = object$info$call,
                    mcall = mcall,
                    formula = formList,
                    direct = object$info$direct,
                    fit = object$info$fit,
                    family = family,
                    control = control,
                    info = stopinfo,
                    model = model,
                    node = nodes,
                    grownode = nodes,
                    nstep = nstep,
                    splitpath = splitpath,
                    dotargs = object$info$dotargs))
    class(tree) <- c("tvcm", "party")
  }

  return(tree)
}


##' -------------------------------------------------------- #
##' \code{\link{tvcm_grow_fit}} fits the current node model.
##'
##' @param call    an object of class call
##' @param doFit   a logical indicating whether the parameters
##'    have to be optimized.
##'
##' @return A list of formulas ('root', 'tree' and 'original').
##'
##' @details Used in 'tvcm' and 'tvcm_grow_sctest'. 'glm.doNotFit'
##' is just an utility function to skip \code{\link{glm.fit}}
##' if \code{doFit = FALSE}
##' -------------------------------------------------------- #

glm.doNotFit <- function(x, y, weights = NULL, start = NULL, etastart = NULL,
                         mustart = NULL, offset = NULL, family = gaussian(),
                         control = list(), intercept = TRUE) {
  coefficients <- rep.int(0, NCOL(x))
  names(coefficients) <- colnames(x)
  if (is.null(weights)) weights <- rep.int(1.0, NROW(x))
  if (is.null(offset)) offset <- rep.int(0.0, NROW(x))
  if (!is.null(start)) {
    if (!is.null(names(start))) {
      start <- start[intersect(names(coefficients), names(start))]
      coefficients[names(start)] <- start
    } else {
      if (length(start) > length(coefficients))
        start <- start[seq_along(coefficients)]
      coefficients[seq_along(start)] <- start
    }
  } 
  return(list(coefficients = coefficients, residuals = NULL,
              effect = NULL, R = NULL, rank = NULL,
              qr = NULL, family = family,
              linear.predictor = etastart,
              deviance = NULL, aic = NULL, null.deviance = NULL,
              iter = 0, weights = NULL, prior.weights = weights,
              df.residual = NULL, df.null = NULL, y = y,
              converged = TRUE, boundary = TRUE))
}

tvcm_grow_fit <- function(mcall, doFit = TRUE) {
    
  ## extract information from 'mcall'
  env <- environment(mcall)

  ## set mcall if coefficients are not to optimized
  if (!doFit) {
    if (inherits(eval(mcall$family, env), "family.olmm")) {
      mcall$doFit <- FALSE
    } else {        
      mcall$method <- glm.doNotFit # skips glm.fit
    }
  }
  
  ## fit model
  object <- suppressWarnings(eval(mcall, env))

  ## return error if fitting failed
  if (inherits(object, "try-error")) stop("model fitting failed.")
  if (doFit && !object$conv) stop("no convergence")
    
  ## return model
  return(object)
}


##' -------------------------------------------------------- #
##' Updates the model matrix and re-fits the current node
##' model. Used for the grid-search in 'tvcm_grow_loss'.
##'
##' @param object a prototype model
##' @param mcall   the mcall for the prototype model
##'
##' @return A list of formulas ('root', 'tree' and 'original').
##'
##' @details Used in 'tvcm' and 'tvcm_grow_loss'. Note that the
##' function will modify the slots of the original object as well!
##'
##' To do:
##' Improve performance for non 'olmm' objects
##' -------------------------------------------------------- #

tvcm_grow_update <- function(object) {
  
  if (inherits(object, "olmm")) {
  
    ## set new partition
    data <- model.frame(object)
    nodeVars <- grep("Node[A-Z]", colnames(data), value = TRUE)
    object$frame[, nodeVars] <- data[, nodeVars]

    ## get terms
    termsFeCe <- terms(object, "fe-ce")
    termsFeGe <- terms(object, "fe-ge")

    ## get constrasts
    contrasts <- object$contrasts
    conCe <- contrasts[intersect(names(contrasts), all.vars(termsFeCe))]
    if (length(conCe) == 0L) conCe <- NULL
    conGe <- contrasts[intersect(names(contrasts), all.vars(termsFeGe))]    
    if (length(conGe) == 0L) conGe <- NULL
    
    ## update model matrix  
    object$X <-
      olmm_merge_mm(x = model.matrix(termsFeCe, object$frame, conCe),
                    y = model.matrix(termsFeGe, object$frame, conGe), TRUE)

    ## extract interaction predictors to be centered
    subsCe <- which(rownames(attr(termsFeCe, "factors")) %in% c("Left", "Right"))
    if (any(subsCe)) {
      subsCe <-
        which(colSums(attr(termsFeCe, "factors")[subsCe,,drop = FALSE]) > 0 &
              !colnames(attr(termsFeCe, "factors")) %in% c("Left", "Right"))
      subsCe <-
        which(attr(object$X, "assign") %in% subsCe & attr(object$X, "merge") == 1)
    }
    subsGe <- which(rownames(attr(termsFeGe, "factors")) %in% c("Left", "Right"))
    if (any(subsGe)) {
      subsGe <-
        which(colSums(attr(termsFeGe, "factors")[subsGe,,drop = FALSE]) > 0 &
              !colnames(attr(termsFeGe, "factors")) %in% c("Left", "Right"))
      subsGe <-
        which(attr(object$X, "assign") %in% subsGe & attr(object$X, "merge") == 2)
    }
    
    ## center the predictors
    for (v in c(subsCe, subsGe)) object$X[,v] <- object$X[,v] - mean(object$X[,v])
    
    ## prepare optimization
    optim <- object$optim
    optim[[1L]] <- object$coefficients
    optim[[4L]] <- object$restricted
    environment(optim[[2L]]) <- environment()
    if (!object$dims["numGrad"]) environment(optim[[3L]]) <- environment()
    FUN <- optim$fit
    subs <- which(names(optim) == "fit")
    optim <- optim[-subs]

    ## run optimization
    output <- try(suppressWarnings(do.call(FUN, optim)), TRUE)
    
    ## check optimized model 
    if (!inherits(output, "try-error")) {
      object$output <- output
      object$conv <- switch(object$optim$fit,
                            optim = object$output$convergence == 0,
                            nlminb = object$output$convergence == 0,
                            ucminf = object$output$convergence %in% c(1, 2, 4))
      if (!object$conv) object <- try(stop("not converged"), TRUE)
    } else {
      object <- output
    }
    
  } else {

    ## modify components in 'object'
    x <- model.matrix(object$formula, model.frame(object))
    object <- try(suppressWarnings(
                    glm.fit(x = x, y = object$y, weights = object$prior.weights,
                            start = object$coefficients, offset = object$offset,
                            family = object$family, control = object$control,
                            intercept = TRUE)), TRUE)
    if (!inherits(object, "try-error")) {
      class(object) <- c("glm", "lm")
      if (!object$conv) object <- try(stop("not converged"), TRUE)
    }
  }
  
  ## return model
  return(object)
}


tvcm_grow_gefp <- gefp.olmm # see 'olmm-methods.R'


##'-------------------------------------------------------- #
##' Computes candidate splits for the current step
##'
##' @param splits   a list. The former list of splits
##' @param partid   a vector of candidate partitions for splitting.
##' @param spart    integer scalar. The partition in which the
##'                 last split was employed
##' @param varid    a \code{list} with a vector for each partition that
##'    that specifies candidate variables for splitting.
##' @param nodeid   a \code{list} with a vector for each partition that
##'    that specifies candidate nodes for splitting.
##' @param model    a fitted model of class \code{\link{olmm}} or
##'    \code{\link{glm}}.
##' @param nodes    a \code{list} with a \code{\link{partynode}}
##'    object for each partition.
##' @param where    a \code{list} with a factor vector for each
##'    partition that the assigns observations to nodes.
##' @param partData a data frame with variables for
##'    splitting.
##' @param control    a \code{list} of control parameters as produced
##'    by 'tvcm_control.'
##' 
##' @return A list of splits. Entries for splits that
##'    exceed the tuning parameters are a vector of length
##'    zero.
##'
##' @details Used in 'tvcm'.
##'-------------------------------------------------------- #

tvcm_grow_setsplits <- function(splits, spart, partid,
                                nodeid, varid, model, nodes,
                                where, partData, control) {
  
  ## get tree size criteria of current tree(s)
  width <- sapply(nodes, width)  
  depth <- lapply(nodes, function(node) {
    unlist(nodeapply(node, nodeids(node, terminal = TRUE), function(node) {
      info_node(node)$depth }))
  })
  w <- weights(model)

  ##'------------------------------------------------------ #
  ##' 'getSplits' extracts splits for each combination of
  ##' partition, variable and node.
  ##'
  ##' @param pid integer. The partition identification number.
  ##' @param vid integer. The row of the variable in 'partData'.
  ##' @param nid integer. The node identification number.
  ##'    For example, if the tree has terminal nodes with
  ##'    identification numbers 2, 4 and 5, the index 1 is used
  ##'    to get splits in node 2, index 2 for node 4 and so on.
  ##'
  ##' @return a matrix with one row for each split.
  ##'------------------------------------------------------ #
  
  getSplits <- function(pid, nid, vid) {

    ## get subset
    subs <- where[[pid]] == levels(where[[pid]])[nid]

    ## get partitioning variable
    z <- partData[, vid]
    
    ## return 'NULL' if tree classical tree growing parameters exceeded
    if (width[pid] >= control$maxwidth[pid] |
        depth[[pid]][nid] >= control$maxdepth[pid] |
        sum(subs) < 1L |
        sum(w[subs]) < 2 * control$minsize[pid]) {
      rval <- matrix(, 0, ifelse(is.numeric(z), 3L, nlevels(z) + 2L))
      colnames(rval) <- c(if (is.numeric(z)) "cut" else levels(z),
                          "loss", "df")
      attr(rval, "type") <- "loss"
      attr(rval, "keeplosscount") <- 0
      return(rval)
    }
    type <- "loss"      

    if (is.numeric(z)) {
      ## continuous variables
      
      sz <- z[subs]
      
      ## get all splits that satisfy the 'minsize' criteria
      subsL <- rev(which(cumsum(w[subs][order(sz)]) < control$minsize[pid]))[1L]
      subsR <- which(cumsum(w[subs][order(sz)]) >
                     (sum(w[subs]) - control$minsize[pid]))[1L]
      
      if (subsL <= subsR && sort(sz)[subsL] < sort(sz)[subsR]) {
        sz <- sort(sz)
        sz <- sz[sz >= sz[subsL] & sz < sz[subsR]]
        rval <- unique(sz)
        
        ## reduce the number of splits according to 'control$maxevalsplit'
        if ((length(rval) - 1)  > control$maxnumsplit) {
          nq <- control$maxnumsplit - 1L
          rval <- c()
          while (length(rval) < (control$maxnumsplit + 1L)) {
            nq <- nq + 1L
            rval <- unique(quantile(sz, c((1:nq) / (nq + 1L)), 1), type = 1L)
          }
          rval <- rval[-length(rval)]
        }
        rval <- cbind(rval, rep.int(NA, length(rval)), rep.int(NA, length(rval)))
        
      } else {
        rval <- matrix(, 0L, 3L)
      }
      colnames(rval) <- c("cut", "loss", "df")
      
    } else if (is.factor(z)) { # categorical variables
      
      nl <- nlevels(z) # all levels
      nld <- nlevels(droplevels(z[subs])) # observed levels in current node
      zdlev <- which(levels(z) %in% levels(droplevels(z[subs])))
      
      if (is.ordered(z)) {
        ## ordinal variables
        
        rval <- diag(nl)
        rval[lower.tri(rval)] <- 1L
        rval <- rval[-nl,, drop = FALSE]
        
      } else {        
        ## nominal variables
        
        if (nld <= control$maxfacsplit) {
          
          ## exhaustive search         
          mi <- 2^(nld - 1L) - 1L
          rval <- .Call("tvcm_nomsplits",
                        as.integer(nl),
                        as.integer(zdlev), as.integer(nld),
                        as.integer(mi), PACKAGE = "vcrpart")
          rval <- matrix(rval, mi, nl, byrow = TRUE)
          
        } else  {
          
          ## Heuristic reduction of splits: in tvcm_grow_loss,
          ## the 'isolated' coefficients of each category are
          ## computed. The coefficients are used for ordering
          ## the categories and finally the variable is treated
          ## as ordinal. See tvcm_grow_loss          
          rval <- diag(nl)
          type <- "coef"
        }
      }
      
      ## delete indadmissible splits
      valid <- tvcm_setsplits_validcats(rval, z, w, subs, control$minsize[pid])
      rval <- rval[valid,, drop = FALSE]
      
      ## delete ordinal splits if 'maxordsplit' is exceeded
      if (is.ordered(z) && nrow(rval) > control$maxordsplit) {         
        nCat <- apply(rval, 1, function(x) sum(subs & z %in% levels(z)[x > 0L]))
        nCat <- round(c(nCat[1L], diff(nCat)))
        zd <- rep.int(1:nrow(rval), nCat)
        nq <- control$maxordsplit - 1L; rows <- 1;
        while(length(rows) < control$maxordsplit) {
          nq <- nq + 1L
          rows <- unique(quantile(zd, (1:nq) / (nq + 1L)), type = 1L)
        }
        rval <- rval[rows,,drop=FALSE]
      }
      
      rval <- cbind(rval, rep.int(NA, nrow(rval)), rep.int(NA, nrow(rval)))
      colnames(rval) <- c(levels(z), "loss", "df")
      
    } else {
      
      rval <- matrix(, 0L, 3L)
      colnames(rval) <- c("cut", "loss", "df")
    }
    attr(rval, "type") <- type
    attr(rval, "keeplosscount") <- 0
    return(rval)
  }

  ## compute splits in all partitions, partitioning variables and nodes  
  rval <- vector("list", length(partid))
  for (pid in seq_along(partid)) {
    rval[[pid]] <- vector("list", length(nodeid[[partid[pid]]]))
    for (nid in seq_along(nodeid[[partid[pid]]])) {
      rval[[pid]][[nid]] <- vector("list", length(varid[[partid[pid]]]))
      for (vid in seq_along(varid[[partid[pid]]])) {
        split <- splits[[pid]][[nid]][[vid]]
        if (is.null(split) | !is.null(split) && attr(split, "type") == "coef" |
            width[pid] >= control$maxwidth[pid]) {
          split <- getSplits(partid[pid],
                             nodeid[[partid[pid]]][nid],
                             varid[[partid[pid]]][vid])
        } else {
          if (nrow(split) > 0L &&
              (spart != partid[pid] |
               attr(split, "keeplosscount") >= control$keeploss)) {
            attr(split, "keeplosscount") <- 0
            split[, c("loss", "df")] <- NA
          } else {
            attr(split, "keeplosscount") <- attr(split, "keeplosscount") + 1L
          } 
        }
        rval[[pid]][[nid]][[vid]] <- split
      }
    }
  }

  ## return list of updated 'splits'
  return(rval)
}


##'------------------------------------------------------ #
##' Computes the valid categorical splits defined in a
##' integer matrix with respect to the \code{minsize}
##' criteria 
##'
##' @param cp an integer matrix that defines the splits
##'    according to the levels in \code{z}.
##' @param z the categorical (nominal or ordered) partitioning
##'    variable.
##' @param weights a weights vector 
##' @param subs a logical vector that extracts the current
##'    node.
##' @param minsize the minimum number of splits in a node
##'
##' @return a logical vector of length equal to the number
##'    of rows of \code{cp}.
##'------------------------------------------------------ #

tvcm_setsplits_validcats <-  function(cp, z, weights, subs, minsize) {

  ## delete cuts that do not satisfy 'minsize'
  rval <- rep.int(TRUE, nrow(cp))
  for (i in seq_along(rval)) {
    Node <- factor(1 * (z[subs] %in% levels(z)[cp[i,]==1L]))
    if (nlevels(Node) == 1L | (nlevels(Node) > 1L &&
                 any(tapply(weights[subs], Node, sum) < minsize)))
      rval[i] <- FALSE
  }

  ## return logical vector
  return(rval)
}


##'------------------------------------------------------ #
##' Updates the list of splits after coefficient
##' constancy tests.
##'
##' @param splits  a 'list' as produced by 'tvcm_grow_splits'
##' @param partid     a vector of candidate partitions for splitting.
##' @param spart      the selected partition
##' @param nodeid     a list with a vector for each partition that
##'    that specifies candidate nodes for splitting.
##' @param snode      the selected node
##' @param varid      a list with a vector for each partition that
##'    that specifies candidate variables for splitting.
##' @param svar       the selected variable

##' 
##' @return An modified list of splits.
##'
##' @details Used in 'tvcm'.
##'------------------------------------------------------ #

tvcm_setsplits_sctest <- function(splits, partid, spart,
                                  nodeid, snode, varid, svar) {

  ## set loss of not selected parts to -Inf
  for (pid in seq_along(partid))
    for (nid in seq_along(nodeid[[pid]]))
      for (vid in seq_along(varid[[pid]])) {
        if (pid == spart & nid == snode & vid == svar) {
          splits[[pid]][[nid]][[vid]][, "loss"] <- NA
        } else {
          splits[[pid]][[nid]][[vid]][, "loss"] <- -Inf
        }
      }
  ## return updated 'splits'
  return(splits)
}


##'------------------------------------------------------ #
##' Updates the list of splits after grid search
##'
##' @param splits  a 'list' as produced by 'tvcm_grow_splits'
##' @param partid     a vector of candidate partitions for splitting.
##' @param spart      the selected partition
##' @param nodeid     a list with a vector for each partition that
##'    that specifies candidate nodes for splitting.
##' @param snode      the selected node
##' @param varid      a list with a vector for each partition that
##'    that specifies candidate variables for splitting.
##' @param svar       the selected variable
##' 
##' @return An modified list of splits.
##'
##' @details Used in 'tvcm'.
##'
##' To do:
##' 2014-07-21: find a better rule! (function is not used currently!)
##'------------------------------------------------------ #

tvcm_setsplits_splitnode <- function(splits, spart, snode,
                                     nodeid, loss, model, control) {
  
  ## expand the splits list
  lnodes <- nodeid[[spart]][nodeid[[spart]] < snode]
  unodes <- nodeid[[spart]][nodeid[[spart]] > snode]
  split0 <- splits[[spart]]
  split <- vector("list", length(split0) + 1L) 
  if (length(lnodes) > 0L) split[lnodes] <- split0[lnodes]
  if (length(unodes) > 0L) split[unodes + 1L] <- split0[unodes]
  splits[[spart]] <- split
  
  ## return updated 'splits'
  return(splits)
}


##'------------------------------------------------------ #
##' Randomly select partitions, variables and nodes
##'
##' @param splits  a 'list' as produced by 'tvcm_grow_splits'
##' @param partid     a vector of candidate partitions for splitting.
##' @param spart      the selected partition
##' @param varid      a list with a vector for each partition that
##'    that specifies candidate variables for splitting.
##' @param svar       the selected variable
##' @param nodeid     a list with a vector for each partition that
##'    that specifies candidate nodes for splitting.
##' @param snode      the selected node
##' 
##' @return An modified list of splits.
##'
##' @details Used in 'tvcm'.
##'------------------------------------------------------ #

tvcm_setsplits_rselect <- function(splits, partid, nodeid, varid, control) {

  ## get the node partitions
  nodeidC <- nodeid
  for (pid in seq_along(partid))
    for (nid in seq_along(nodeid[[pid]]))
      if (all(sapply(splits[[pid]][[nid]], function(x) length(x) == 0L)))
        nodeidC[[pid]] <- setdiff(nodeidC[[pid]], nid)

  ## get the partition candidates
  partidC <- partid
  for (pid in seq_along(partid))
    if (length(nodeidC[[pid]]) == 0L) setdiff(partidC, pid)

  ## get variable candidates for each partition
  varidC <- varid
  for (pid in seq_along(partid))
    for (vid in seq_along(varid[[pid]]))
      if (all(sapply(splits[[pid]], length) == 0L))
        varidC[[pid]] <- setdiff(varidC[[pid]], vid)


  ## random selections
  spart <- sort(sample(partidC, min(length(partidC), control$ptry)))
  svar <- lapply(seq_along(varidC), function(pid) {
    s <- sample(length(varidC[[pid]]), min(length(varidC[[pid]]), control$vtry[pid]))
    return(varidC[[pid]][sort(s)])
  })
  snode <- lapply(seq_along(nodeidC), function(pid) {
    s <- sample(length(nodeidC[[pid]]), min(length(nodeidC[[pid]]),control$vtry[pid]))
    return(nodeidC[[pid]][sort(s)])
  })

  ## delete not selected nodes from 'splits'
  for (pid in seq_along(partid)) 
    for (nid in seq_along(nodeid[[pid]])) 
      for (vid in seq_along(varid[[pid]])) 
        if (nrow(splits[[pid]][[nid]][[vid]]) > 0 &&
            !(pid %in% spart & vid %in% svar[[pid]] & nid %in% snode[[pid]]))
          splits[[pid]][[nid]][[vid]][, "loss"] <- -Inf

  ## return updated 'splits'
  return(splits)
}


##'------------------------------------------------------ #
##' Processing of nodewise coefficient constancy tests.
##'
##' @param model    the current model.
##' @param nodes    an object of class 'partynode'.
##' @param where    a list of vectors that locate the observations
##'    to their corresponding node(s)
##' @param partid   an integer vector that indicates which
##'    partitions should be tested.
##' @param varid    a list with a vector for each partition that
##'    that specifies the variables to be tested
##' @param partData a 'data.frame' with the partitioning variables.
##' @param control  an object of class 'tvcm_control'.
##' 
##' @return A list with partitions 'statistic' and 'p.value'.
##'
##' @details Used in 'tvcm'.
##'------------------------------------------------------ #

tvcm_grow_sctest <- function(model, nodes, where, partid, nodeid, varid, 
                            splits, partData, control) {
    
  ## get variable types
  functional <- sapply(partData, function(x) {
    switch(class(x)[1],
           factor = control$functional.factor,
           ordered = control$functional.ordered,
           integer = control$functional.numeric,
           numeric = control$functional.numeric)
  })
  
  ## prepare list with arguments for 'sctest'
  rval <- vector("list", length(nodes))
  for (pid in seq_along(nodes)) {
    dim <- c(nlevels(where[[pid]]), length(varid[[pid]]), control$ninpute)
    dn <- list(paste("Node", LETTERS[pid], levels(where[[pid]]), sep = ""),
               colnames(partData)[varid[[pid]]], 1:control$ninpute)
    rval[[pid]] <- array(, dim = dim, dimnames = dn)
  }

  ## call 'estfun'
  eCall <- list(name = as.name(ifelse(inherits(model, "olmm"),"estfun.olmm", "estfun")))
  eCall$x <- quote(model)
  eCall[names(control$estfun)] <- control$estfun
  mode(eCall) <- "call"
  scores <- replicate(control$ninpute, eval(eCall))

  ## set the 'gefp' call (which is called in each iteration below)
  gCall <- call(name = "tvcm_grow_gefp", object = quote(model),
                scores = quote(sc),
                order.by = quote(z), subs = quote(rows),
                parm = quote(cols), center = TRUE, silent = TRUE)

  terms <- # useful information to identify the coefficients to test
    tvcm_get_terms(dimnames(scores)[[2L]],
                   lapply(nodes, function(node) nodeids(node, terminal = TRUE)),
                   control$parm)
  
  ## apply test for each variable and partition separately
  for (pid in seq_along(partid)) { # loop over partitions
      
    for (vid in seq_along(varid[[pid]])) { # loop over variables
      
      ## get variable to test
      z <- partData[, varid[[pid]][vid]]       
      
      for (nid in seq_along(nodeid[[pid]])) { # loop over nodes
        
        ## check if there is a permitted split
        if (length(splits[[pid]][[nid]][[vid]]) > 0L) {
          
          ## observations of the current partition
          rows <- where[[pid]] == levels(where[[pid]])[nid] 
          
          ## columns corresponding to the tested subset
          if (nlevels(where[[pid]]) == 1L) {
            cols <- unlist(control$parm[[pid]])
          } else {
            cols <- terms$partition == LETTERS[pid] & terms$node == 
              levels(where[[pid]])[nid]
            cols <- dimnames(scores)[[2L]][cols]
          }
          
          for (k in 1:control$ninpute) {
            sc <- matrix(scores[,,k,drop=FALSE], dim(scores)[1], dim(scores)[2],
                         dimnames = dimnames(scores)[1L:2L])            
            gefp <- try(eval(gCall), TRUE)
            
            ## extract test statistic
            if (!inherits(gefp, "try-error")) {            
              if (is.character(functional)) {
                functional <- tolower(functional)
                fi <- switch(functional[varid[[pid]][vid]],
                             "suplm" = supLM(from = control$trim),
                             "lmuo" = catL2BB(gefp),
                             stop("Unknown efp functional."))
              } else {
                fi <- functional[varid[[pid]][vid]]
              }              
              test <- try(sctest(x = gefp, functional = fi), TRUE)
              if (!inherits(test, "try-error")) {
                rval[[pid]][nid, vid, k] <- test$p.value
              }
            }
          }
        }
      }
    }
  }
  rval <- lapply(rval, function(x) apply(x, c(1L, 2L), mean, 
                                         na.rm = TRUE))
  rval <- lapply(rval, function(x) {
      x[is.nan(x)] <- NA
      return(x)
  })
  names(rval) <- LETTERS[partid]
  return(rval)
}


tvcm_sctest_bonf <- function(test, type) {
  for (pid in seq_along(test)) {
    for (nid in 1:nrow(test[[pid]])) {
      k <- switch(type,
                  "none" = 1.0,
                  "all" = sum(!is.na(unlist(test))),
                  "nodewise" = sum(!is.na(test[[pid]][nid, ])),
                  "partitionwise" = sum(!is.na(test[[pid]])))
      pval1 <- pmin(1.0, k * test[[pid]][nid, ])
      pval2 <- c(1.0 - (1.0 - test[[pid]][nid, ]) ^ k)
      test[[pid]][nid, ] <-
        ifelse(!is.na(test[[pid]][nid, ]) & (test[[pid]][nid, ] > 0.01),
               pval2, pval1)
    }
  }
  return(test)
}


##'-------------------------------------------------------- #
##' Computes the loss for each possible split.
##'
##' @param varid    an integer vector indicating the partitioning
##'    variables to be evaluated.
##' @param nodeid   an integer vector indicating the nodes to be
##'    evaluated.
##' @param model    the current model
##' @param nodes    an object of class 'partynode'.
##' @param partData a 'data.frame' with the partitioning variables.
##' @param control  an object of class 'tvcm_control'.
##' @param mcall
##' @param step     integer. The current step number.
##'
##' @return A nested list with loss matrices. Partitions of nodes
##'    are nested in partitions for variables. 
##'
##' @details Used in 'tvcm'.
##'-------------------------------------------------------- #

tvcm_grow_loss <- function(splits, partid, nodeid, varid, 
                           model, nodes, where, partData, 
                           control, mcall, formList, step) {

  verbose <- control$verbose; control$verbose <- FALSE;

  loss0 <- control$lossfun(model)
  mfName <- switch(deparse(mcall[[1]]), glm = "model", olmm = "frame")
  
  if (verbose) cat("\n* computing splits ")
  
  ## function to get the splitstatistic for a given cutpoint
  getSplitstat <- function(cutpoint, type = "loss",
                           pid, nid, vid, 
                           model, modelNuis,
                           nuisance) {
    
    ## set node indicator
    subs <- where[[pid]] == levels(where[[pid]])[nid]
    z <- partData[, vid]    
    if (is.numeric(z)) {
      zs <- z <= cutpoint
    } else {
      zs <- z %in% levels(z)[cutpoint > 0L]            
    }
    model[[mfName]]$Left <- 1 * (subs & zs)
    model[[mfName]]$Right <- 1 * (subs & !zs)

    ## fit the 'update' model
    model <- tvcm_grow_update(model)
    
    if (type == "loss") {
      if (!inherits(model, "try-error")) {
        rval <- c(loss0 - control$lossfun(model),
                  length(coef(model)[grep("Left", names(coef(model)))]) -
                  length(nuisance))
        if (is.null(modelNuis)) {
          return(rval)
        } else {
          modelNuis[[mfName]]$Left <- 1 * (subs & zs)
          modelNuis[[mfName]]$Right <- 1 * (subs & !zs)
          modelNuis <- tvcm_grow_update(modelNuis)
          rval[1L] <- rval[1L] - (loss0 - control$lossfun(modelNuis))
          return(rval)
        }  
      } else {
        return(c(NA, NA))        
      }
    } else {
      parm <- grep("Left", names(coef(model)), value = TRUE)
      nuis <- gsub("Node[A-Z]", "Left", nuisance)
      parm <- setdiff(parm, nuis)  
      if (inherits(model, "try-error")) {
        return(rep.int(NA, length(parm)))
      } else {
        return(coef(model)[parm])
      }
    }
  }

  mcall$data <- eval(mcall$data, environment(mcall))
  w <- weights(model)
  mcall$offset <- predict(model, type = "link")
  if (inherits(model, "glm")) {
    mcall$x <- TRUE
    mcall$y <- TRUE
    mcall$model <- TRUE
  } else if (inherits(model, "olmm")) {
    mcall$restricted <- grep("ranefCholFac", names(coef(model)), value = TRUE)
    mcall$start <- coef(model)[mcall$restricted]
  }

  ff <- tvcm_formula(formList, rep.int(FALSE, length(partid)), 
                     eval(mcall$family, environment(mcall)),
                     environment(mcall), full = FALSE, update = TRUE)
  Left <- sample(c(0, 1), nobs(model), replace = TRUE)
  Right <- Left - 1
  
  for (pid in seq_along(partid)) { 

    if (length(unlist(splits[[pid]])) > 0L) {
     
      mcall$formula <- ff$update[[pid]][[1L]]
      mcall$data$Left <- Left
      mcall$data$Right <- Right 
      sModel <- tvcm_grow_fit(mcall, doFit = FALSE)
      sModel$coefficients[grepl("Left", names(sModel$coefficients))] <- 0.0
      sModel$coefficients[grepl("Right", names(sModel$coefficients))] <- 0.0
      sModel$control <- model$control
      
      if (length(control$nuisance[[pid]]) == 0L) {
        sModelN <- NULL
      } else {
        mcallN <- mcall
        mcallN$formula <- ff$update[[pid]][[2L]]
        sModelN <- tvcm_grow_fit(mcallN, doFit = FALSE)
        sModelN$coefficients[] <- 0.0
        sModelN$control <- model$control
      } 
      
      ## run computation
      for (nid in seq_along(splits[[pid]])) {
        if (length(unlist(splits[[pid]][[nid]])) > 0L) {
          for (vid in seq_along(splits[[pid]][[nid]])) {
            cp <- splits[[pid]][[nid]][[vid]]
            type <- attr(cp, "type")
            subs <- is.na(cp[, "loss"])
            cp <- cp[, !colnames(cp) %in% c("loss", "df"), drop = FALSE]
            if (any(subs)) {
              st <- apply(cp, 1, getSplitstat, type = type,
                          pid = partid[pid],
                          nid = nodeid[[partid[pid]]][nid],
                          vid = varid[[partid[pid]]][vid],
                          model = sModel, modelNuis = sModelN,
                          nuisance = control$nuisance[[pid]])
              if (is.matrix(st)) st <- t(st) else st <- matrix(st, ncol = 1L)

              if (type == "loss") {
                splits[[pid]][[nid]][[vid]][subs, c("loss", "df")] <- st
                
              } else if (type == "coef") {
                
                ## if 'z' is a nominal variable with many categories               
                z <- partData[, varid[[partid[pid]]][vid]]
                subs <- where[[partid[pid]]] ==
                  levels(where[[partid[pid]]])[nodeid[[partid[pid]]][nid]]
                
                invalids <- attr(na.omit(st), "na.action")
                if (length(invalids) > 0L) {
                  cp <- cp[-invalids,,drop = FALSE]
                  st <- st[-invalids,,drop = FALSE]
                }
                score <- rep.int(0, nlevels(z))
                score[colSums(cp) > 0] <- prcomp(st)$x[,1]
                
                ## define 'z' as ordinal and retrieve the splits             
                zd <- factor(z, levels = levels(z)[order(score)], ordered = TRUE) 
                nl <- nlevels(zd)
                cp <- diag(nl)
                cp[lower.tri(cp)] <- 1L
                cp <- cp[-nl,, drop = FALSE]
                valid <-
                  tvcm_setsplits_validcats(cp, zd, w, subs,
                                           control$minsize[partid[pid]])
                cp <- cp[valid,, drop = FALSE]
                
                ## reorder the columns of 'cp' acc. to the original categories
                cp <- cp[,order(order(score)), drop = FALSE]
                
                if (nrow(cp) > control$maxordsplit) {
                  nCat <-
                    apply(cp, 1, function(x) sum(subs & z %in% levels(z)[x > 0L]))
                  nCat <- round(c(nCat[1L], diff(nCat)))
                  zd <- rep.int(1:nrow(cp), nCat)
                  nq <- control$maxordsplit - 1L; rows <- 1;
                  while(length(rows) < control$maxordsplit) {
                    nq <- nq + 1L
                    rows <- unique(quantile(zd, (1:nq) / (nq + 1L)), type = 1L)
                  }
                  cp <- cp[rows,,drop=FALSE]
                }
                
                ## compute the loss of the new splits
                st <- apply(cp, 1, getSplitstat, type = "loss",
                            pid = partid[pid],
                            nid = nodeid[[partid[pid]]][nid],
                            vid = varid[[partid[pid]]][vid],
                            model = sModel,
                            modelNuis = sModelN,
                            nuisance = control$nuisance[[pid]])

                if (is.matrix(st)) st <- t(st) else st <- matrix(st, ncol = 2L)
                split <- cbind(cp, st)
                colnames(split) <- c(levels(z), "loss", "df")
                attr(split, "type") <- "coef"
                splits[[pid]][[nid]][[vid]] <- split
              }
            }
          }
        }
      }
    }
  }
  
  ## function that extracts the loss reduction (eventually corrected by the
  ## number of predictors)
  getLossDiff <- function(x) {
    if (is.list(x)) return(lapply(x, getLossDiff))
    if (is.matrix(x))
        if (nrow(x) > 0L)
            return(x[, "loss"] - control$dfpar * x[, "df"]) else return(numeric())
    return(x)
  }
  loss <- getLossDiff(splits)
  
  ## function that extracts the maximum loss reduction
  getMaxLossDiff <- function(x) {
      x <- unlist(x)
      if (length(x) == 0L) return(-Inf)
      x <- na.omit(x)
      if (length(x) > 0L) return(max(x)) else return(-Inf)
  }

  maxLossDiff <- max(c(-Inf, na.omit(unlist(loss))))
  
  if (maxLossDiff > -Inf) {
    
    ## select the partition, node and variable
    spart <- which(sapply(sapply(loss, getMaxLossDiff), identical, maxLossDiff))
    if (length(spart) > 1L) spart <- sample(spart, 1L)
    snode <-
      which(sapply(sapply(loss[[spart]], getMaxLossDiff), identical, maxLossDiff))
    if (length(snode) > 1L) snode <- sample(snode, 1L)
    svar <- which(sapply(sapply(loss[[spart]][[snode]], getMaxLossDiff),
                         identical, maxLossDiff))
    if (length(svar) > 1L) svar <- sample(svar, 1L)
    
    ## select the cut
    stat <- splits[[spart]][[snode]][[svar]]
    cutid <- which(stat[, "loss"] == max(stat[, "loss"]))
    if (length(cutid) > 1L) cutid <- sample(cutid, 1L)
    
    if (verbose) cat("OK")
    
    return(list(partid = partid[spart],
                nodeid = nodeid[[partid[spart]]][snode],
                varid = varid[[partid[spart]]][svar],
                cutid = cutid,
                cut = stat[cutid, !colnames(stat) %in% c("loss", "df")],
                loss = as.numeric(stat[cutid, "loss"]),
                ploss = maxLossDiff,
                df = as.numeric(stat[cutid, "df"]),
                lossgrid = splits))
  } else {
    
    if (verbose) cat("failed")
    
    return(list(partid = NULL, nodeid = NULL, varid = NULL, 
                cutid = NULL, cut = NULL, loss = NULL,
                lossgrid = splits))
    
  }
}


##'-------------------------------------------------------- #
##' Incorporates a new binary split into an existing
##' tree structire.
##'
##' @param nodes    an object of class 'partynode'.
##' @param loss     a list produced by 'tvcm_grow_loss'.
##' @param partData a 'data.frame' with the partitioning variables.
##' @param step     integer. The current algorithm step.
##'
##' @return A list of formulas ('root', 'tree' and 'original').
##'
##' Used in 'tvcm'.
##'-------------------------------------------------------- #

tvcm_grow_splitnode <- function(nodes, where, loss, partData, step, weights) {

  pid <- loss$partid
  nid <- loss$nodeid
  vid <- loss$varid
  nidLab <- nodeids(nodes[[pid]], terminal = TRUE)[nid]
  stat <- loss$loss
  cut <- loss$cut
  x <- partData[, vid]
  
  ## collect information for the split
  subs <- where[[pid]] == levels(where[[pid]])[nid]
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
  oldnodes <- as.list(nodes[[pid]])
  
  ## setup 'newnodes' object
  subsN <- which(sapply(oldnodes, function(node) node$id) == nidLab)
  newnodes <- vector("list", length(oldnodes) + 2)
  newnodes[1:subsN] <- oldnodes[1:subsN]
  if (length(oldnodes) > nidLab)
    newnodes[(subsN + 3L):length(newnodes)] <-
      oldnodes[(subsN + 1L):length(oldnodes)]
  
  ## adjust ids of children
  ids <- sapply(newnodes, function(x) if (!is.null(x$id)) x$id else -1)
  for (i in 1L:length(newnodes))
    if (!is.null(newnodes[[i]]$kids))
      newnodes[[i]]$kids <- which(ids %in% newnodes[[i]]$kids)
  
  ## setup new split
  newnodes[[subsN]]$split <-
    partysplit(varid = vid, breaks = breaks, index = index,
               info = list(ordered = ordered, step = step))
  newnodes[[subsN]]$kids <- nidLab + 1L:2L
  newnodes[[subsN]]$info$dims <- c(n = sum(weights[subs]))
  
  ## add new children
  newnodes[[subsN + 1L]] <-
    list(id = nidLab + 1L,
         info = list(
           dims = c(n = sum(weights[subsL])),
           depth = newnodes[[subsN]]$info$depth + 1L))
  
  newnodes[[subsN + 2L]] <-
    list(id = nidLab + 2L,
         info =
         list(dims = c(n = sum(weights[subsR])),
              depth = newnodes[[subsN]]$info$depth + 1L))
  
  ## adjust ids
  for (i in 1L:length(newnodes))
    newnodes[[i]]$id <- i
  
  ## return new nodes
  nodes[[pid]] <- as.partynode(newnodes)
  return(nodes)
}


##'-------------------------------------------------------- #
##' Extracts the formula for the root node and the tree
##' from the output of \code{\link{vcrpart_formula}}.
##'
##' @param formList a list of formulas from 'vcrpart_formula'.
##' @param root     logical vector of the same length as the
##'    'vc' slot of 'formList'.
##' @param family   an object of class 'family' or 'family.olmm'.
##' @param env      the environment where the output formula
##'    is to be evaluated.
##' @return A list of formulas ('root', 'tree' and 'original').
##'
##' @details Used in \code{\link{predict.fvcm}} and
##'    \code{\link{tvcm}}.
##'-------------------------------------------------------- #

tvcm_formula <- function(formList, root, family = cumulative(),
                         env = parent.frame(),
                         full = TRUE, update = FALSE) {

  yName <- rownames(attr(terms(formList$original), "factors"))[1L]

  ## puts the predictors for fixed effects and varying effects
  ## into one formula
  getTerms <- function(x, effect, root, family, nuisance = NULL) {
    
    ## get 'vc' terms
    vcTerms <- x$vc
    if (!is.null(vcTerms)) {
      vcTerms <-
        lapply(x$vc, function(x) attr(terms(x$eta[[effect]],
                                            keep.order = TRUE), "term.labels"))
      for (i in 1:length(root)) {
        if (root[i]) {
          vcTerms[[i]] <-
            vcTerms[[i]][vcTerms[[i]] != paste("Node", LETTERS[i], sep = "")]
          vcTerms[[i]] <- gsub("Node[A-Z]:", "", vcTerms[[i]])    
        }
      }
      vcTerms <- unlist(vcTerms)
    }

    ## get 'fe' terms
    feTerms <- x$fe$eta[[effect]]
    if (!is.null(feTerms))
      feTerms <- attr(terms(feTerms, keep.order = TRUE), "term.labels")
    
    vcTerms <- setdiff(vcTerms, feTerms)
    
    rval <- ""
    if (length(vcTerms) > 0L)
      rval <- paste(rval, paste(vcTerms, collapse = "+"), sep = "")
    if (length(vcTerms) > 0L & length(feTerms) > 0L)
      rval <- paste(rval, "+", sep = "")
    if (length(feTerms) > 0L)
      rval <- paste(rval, paste(feTerms, collapse = "+"), sep = "")
    if (rval != "" && inherits(family, "family.olmm"))
      rval <- paste(effect, "(", rval, ")", sep = "")

    
    return(c(vcTerms, feTerms))
  }
  
  ## incorporate fixed effects terms
  feCeTerms <- getTerms(formList, "ce", root, family)
  feGeTerms <- getTerms(formList, "ge", root, family) 
  
  ## intercept
  feInt <- formList$fe$intercept
  vcInt <- unlist(lapply(formList$vc, function(x) x$intercept))
  direct <- sapply(formList$vc, function(x) x$direct)
  if (!is.null(vcInt) && any(direct) && root[direct])
    feInt <- "ce"
  
  ## random effects
  if (!is.null(formList$re)) {
    subjectName <- attr(terms(formList$re$cond), "term.labels")      
    getReTerms <- function(effect) {
      rval <- attr(terms(formList$re$eta[[effect]], keep.order = TRUE), "term.labels")
      if (formList$re$intercept == effect) rval <- c("1", rval)
      if (length(rval) == 0L) return(NULL)
      rval <- paste(rval, collapse = "+")
      if (inherits(family, "family.olmm"))
        rval <- paste(effect, "(", rval, ")", sep = "")
      return(rval)
    }
    reForm <- unlist(lapply(c("ce", "ge"), getReTerms))
    reForm <- paste(reForm, collapse = "+")
    if (inherits(family, "family.olmm")) {        
      reForm <- paste("re(", reForm, "|", subjectName,
                      ",intercept='", formList$re$intercept, "')", sep = "")
    } else {        
      if (formList$re$intercept == "none")
        reForm <- paste(reForm, "-1", sep = "")
      reForm <- paste("(", reForm, "|", subjectName, ")", sep = "")
    }        
  } else {
    reForm <- NULL
  }

  getForm <- function(yName, feCeTerms, feGeTerms, feInt, reForm, family, env) {

    fTree <- "" # the return value
    
    feCeForm <- if (length(feCeTerms) == 0L) "" else paste(feCeTerms, collapse = "+")
    if (feCeForm != "" & inherits(family, "family.olmm"))
      feCeForm <- paste("ce(", feCeForm, ")", sep = "")

    feGeForm <- if (length(feGeTerms) == 0L) "" else paste(feGeTerms, collapse = "+")
    if (feGeForm != "" &  inherits(family, "family.olmm"))
      feGeForm <- paste("ge(", feGeForm, ")", sep = "")

    if (feCeForm != "") fTree <- paste(fTree, feCeForm, sep = "")
    if (feCeForm != "" & feGeForm != "") fTree <- paste(fTree, " + ", sep = "")
    if (feGeForm != "") fTree <- paste(fTree, feGeForm, sep = "")
    if (fTree == "") fTree <- "1"

    if (inherits(family, "family.olmm")) {
      fTree <- paste(fTree, ", intercept='", feInt, "'", sep = "")
    } else {
      if (feInt == "none")
        fTree <- paste("-1", fTree, sep = "+")
    }
    
    if (inherits(family, "family.olmm"))
      fTree <- paste("fe(", fTree, ")", sep = "")

    if (!is.null(reForm)) fTree <- paste(fTree, "+", reForm)

    fTree <- paste(yName, "~", fTree)
    
    return(as.formula(fTree, env = env))
  }

  ## full formula
  
  fFull <- NULL
  if (full)
    fFull <- getForm(yName, feCeTerms, feGeTerms, feInt, reForm, family, env)
  
  ## update formulas
  
  fUpdate <- NULL
  if (update) {

    ## get nuisance terms
    nuisance <- lapply(seq_along(formList$vc), function(pid) {
      return(formList$vc[[pid]]$nuisance)
    })
  
    ## update formulas for tvcm_grow_loss
    fUpdate <- vector("list", length(formList$vc))
    for (pid in seq_along(fUpdate)) {
      fUpdate[[pid]] <- vector("list", 2L)
      nLab <- paste("Node", LETTERS[pid], sep = "")
      
      ## full formula
      feCeTmp <- feCeTerms[grep(nLab, feCeTerms)]
      feCeTmp <- c(gsub(nLab, "Left", feCeTmp), gsub(nLab, "Right", feCeTmp))
      feGeTmp <- feGeTerms[grep(nLab, feGeTerms)]
      feGeTmp <- c(gsub(nLab, "Left", feGeTmp), gsub(nLab, "Right", feGeTmp))
      fUpdate[[pid]][[1L]] <-
        getForm(yName,feCeTmp,feGeTmp,"none",reForm, family, env)
      
      ## null formula
      feCeTmp <- feCeTerms[grep(nLab, feCeTerms)]
      feCeTmp <- feCeTmp[grep(nLab, feCeTmp)]
      feCeTmp <- intersect(feCeTmp, nuisance[[pid]])
      feCeTmp <- c(gsub(nLab, "Left", feCeTmp), gsub(nLab, "Right", feCeTmp))
      feGeTmp <- feGeTerms[grep(nLab, feGeTerms)]
      feGeTmp <- feGeTmp[grep(nLab, feGeTmp)]
      feGeTmp <- intersect(feGeTmp, nuisance[[pid]])
      feGeTmp <- c(gsub(nLab, "Left", feGeTmp), gsub(nLab, "Right", feGeTmp))
      fUpdate[[pid]][[2L]] <- getForm(yName,feCeTmp,feGeTmp,"none",reForm,family, env)
    }
  }
  return(list(full = fFull, update = fUpdate))
}


##'-------------------------------------------------------- #
##' Adds a new slot 'parm' to the 'control_tvcm' object
##' for internal purposes.
##'
##' @param control  an object of class 'tvcm_control'.
##' @param model    a root node regression model, e.g., an 'olmm'
##'    or a 'glm' object
##' @param formList a list of formulas from 'vcrpart_formula'.
##' 
##'
##' @return An updated 'tvcm_control' object.
##'
##' @details Used in 'tvcm'.
##'-------------------------------------------------------- #

tvcm_grow_setcontrol <- function(control, model, formList, root, parm.only = TRUE) {

  family <- model$family
  if (is.null(formList$vc)) return(control)

  ## specify the tree size parameters separately for each partition
  
  if (!parm.only) {

    npart <- length(formList$vc)
    
    if (!length(control$maxwidth) %in% c(1L, npart))
      stop("'maxwidth' must be either of length ", 1L, " or ", npart, ".")
    control$maxwidth <- rep_len(control$maxwidth, npart)
    
    if (!is.null(control$minsize) && !length(control$minsize) %in% c(1L, npart))
      stop("'minsize' must be either of length ", 1L, " or ", npart, ".")
    if (!is.null(control$minsize))
      control$minsize <- rep_len(control$minsize, npart)
    
    if (!length(control$maxdepth) %in% c(1L, npart))
      stop("'maxdepth' must be either of length ", 1L, " or ", npart, ".")
    control$maxdepth <- rep_len(control$maxdepth, npart)
    
    ## update 'fvcm' parameters
    
    if (!length(control$ntry) %in% c(1L, npart))
      stop("'ntry' must be either of length ", 1L, " or ", npart, ".")
      control$ntry <- rep_len(control$ntry, npart)
    
    if (!length(control$vtry) %in% c(1L, npart))
      stop("'vtry' must be either of length ", 1L, " or ", npart, ".")
    control$vtry <- rep_len(control$vtry, npart)
  }
  
  ## update the 'parm' and the 'nuisance' slots
  
  ## set 'vcterms' slot
  vcParm <- lapply(formList$vc, function(x) lapply(x$eta, function(x) attr(terms(x), "term.labels")))

  for (pid in seq_along(vcParm)) {
    for (j in names(vcParm[[pid]])) {
      terms <- vcParm[[pid]][[j]]
      if (root[pid]) {
        terms <- terms[terms != paste("Node", LETTERS[pid], sep = "")]
        terms <- gsub("Node[A-Z]:", "", terms)
      }  
      if (length(terms) > 0L) {
        type <- paste("fe", j, sep = "-")
        X <- model.matrix(model, which = type)
        assign <- attr(X, "assign")
        subs <- which(attr(terms(model, type), "term.labels") %in% terms)
        if (length(subs) == 0L) {
          terms <- sapply(terms, function(x) {
            x <- strsplit(x, ":")[[1L]]
            len <- length(x)
            x <- c(if (len > 2) x[1:(len-2)], x[len], x[len-1])
            return(paste(x, collapse = ":"))
          })
          subs <- which(attr(terms(model, type), "term.labels") %in% terms)
        }
        if (length(subs) > 0L) terms <- colnames(X)[assign %in% subs]
      }
      vcParm[[pid]][[j]] <- terms
    }
  }  
  if (inherits(family, "family.olmm")) {
    for (pid in seq_along(vcParm)) {
      if ((len <- length(vcParm[[pid]][[1L]])) > 0L)
        vcParm[[pid]][[1L]] <-
          paste("Eta", rep(1L:model$dims["nEta"], each = len), ":",
                rep(vcParm[[pid]][[1L]], model$dims["nEta"]), sep = "")
    }
  }  
  control$parm <- vcParm
  
  ## set 'intercept' slot (which is always in the first 'vc' term)
  if (control$direct && root[1L]) {
    if (inherits(model, "olmm")) {
      control$parm[[1L]]$ce <- c(grep("Eta[1-9]+:\\(Intercept\\)",
                                      names(coef(model)), value = TRUE),
                                 control$parm[[1L]]$ce)
    } else {
      control$parm[[1L]]$ce <- c("(Intercept)", control$parm[[1L]]$ce)
    }
  }
  
  ## set 'nuisance' slots
  control$nuisance <- lapply(formList$vc, function(x) x$nuisance)
  control$estfun$nuisance <-
    unique(c(control$estfun$nuisance,
             setdiff(names(coef(model)), unlist(control$parm))))
  return(control)
}


##'-------------------------------------------------------- #
##' Extract the node vector from 'newdata' and assign
##' the contrasts.
##'
##' @param object  a fitted 'tvcm' object.
##' @param newdata a data.frame from which the nodes are to be
##'    extracted.
##' @param weights a numeric vector of weights with the same
##'    size than data.
##' @param formList a list of formulas as produced by
##'    \code{vcrpart_formula}.
##'
##' @return A list with values of the variables in 'data'
##'
##' @details Used in tvcm.predict and tvcm.prune.
##'    \code{tvcm_get_fitted} is merely a help function
##'-------------------------------------------------------- #

tvcm_get_fitted <- function(pid, object, newdata, weights, setContrasts) {
  fitted <- fitted_node(object$info$node[[pid]],
                        newdata[,colnames(object$data), drop = FALSE])
  fitted <- factor(fitted, nodeids(object$info$node[[pid]], terminal = TRUE))
  names(fitted) <- rownames(newdata) 
  if (nlevels(fitted) > 1L) {
    if (setContrasts) {
      fitted <- vcrpart_contr.sum(fitted, weights)
    } else {
      contrasts(fitted) <-
        object$contrasts[, paste("Node", LETTERS[pid], sep = "")]
    }
  }
  return(fitted)
}

tvcm_get_node <- function(object, newdata, setContrasts = FALSE, weights,
                          formList = NULL) {
  if (is.null(formList))
    formList <- vcrpart_formula(object$info$formula, object$info$family)
  fitted <- lapply(seq_along(object$info$node), tvcm_get_fitted, object = object,
                   newdata = newdata, weights = weights, setContrasts = setContrasts)
  names(fitted) <- paste("Node", LETTERS[seq_along(object$info$node)], sep = "")
  return(fitted)
}


##'-------------------------------------------------------- #
##' Creates a list which can be used to extract the varying
##' coefficients.
##'
##' @param names character vector. Names of coefficients of the
##'    current model.
##' @param ids   character vector. Names of the current nodes.
##' @param parm  the 'control$parm' slot
##'
##' @return A list with slots
##'    names:     the original coefficient names.
##'    terms:     the names of the terms to which coefficients
##'               belong according to the original formula.
##'    type:      which type of 'fe', 'vc' or 're' the term
##'               belongs to.
##'    node:      the node to which a coefficient belongs to.
##'    partition: the partition to which a coefficients
##'               belongs to.
##'-------------------------------------------------------- #

tvcm_get_terms <- function(names, ids, parm) {
  
  parm <- lapply(parm, function(x) {
    lapply(x, function(x) x[grepl("Node[A-Z]", x)])
    })
  
  if (any(unlist(ids) == 1L)) {
    getNames <- function(x) {
      if (!x %in% unlist(parm)) return(x)
      pid <- which(sapply(parm, function(p) x %in% unlist(p)))
      if (grepl("(Intercept)", x))
        return(gsub("(Intercept)", paste("Node", LETTERS[pid], 1, sep = ""),
                    x, fixed = TRUE))
      if (!grepl("Node", x)) {
        x <- strsplit(x, ":")
        x <- rep(x, length(pid))
        subs <- ifelse(grepl("Eta[1-9]+", x[[1L]][1L]), 2L, 1L)
        for (i in 1:length(x)) {
          x[[i]][subs] <- paste("Node", LETTERS[pid[i]], 1, ":",
                                x[[i]][subs], sep = "")
          x[[i]] <- paste(x[[i]], collapse = ":")
        }
      }
      return(x)
    }
    names <- unlist(lapply(names, getNames))
  }
  
  nodes <- unlist(lapply(seq_along(ids), function(i)
                         paste("Node", names(ids)[i], ids[[i]], sep = "")))
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
        rval <- substr(x[subs], 6, 500)
      return(rval)
    })
  partition <-
    sapply(split, function(x) {
      rval <- ""
      if (any(subs <- x %in% nodes))
        rval <- substr(x[subs], 5, 5)
      return(rval)
    })
  return(list(names = names, 
              terms = terms, type = type,
              node = node, partition = partition))
}


##'-------------------------------------------------------- #
##' Extracts the names of the predictors on 'vc' terms
##' 
##' @param object a 'tvcm' object.
##'
##' @return A character vector.
##'-------------------------------------------------------- #

tvcm_get_vcparm <- function(object) {
    vcTerms <- terms(object$info$formula$original, specials = "vc")
    vcTerms <- rownames(attr(vcTerms, "factors"))[attr(vcTerms, "specials")$vc]
    parm <- lapply(vcTerms, function(x) {
        eta <- eval(parse(text = x))$eta
        if (inherits(object$info$model, "olmm")) {
            etaList <- vcrpart_formula(eta)
            parmCe <- all.vars(etaList$fe$eta$ce)
            parmCe <- paste("Eta", 1:object$info$model$dims["nEta"], ":", parmCe, sep = "")
            parmGe <- all.vars(etaList$fe$eta$ge)
            parm <- c(parmCe, parmGe)
        } else {
            parm <- all.vars(eta)
        }
        return(parm)
    })
    for (pid in seq_along(object$info$formula$vc)) {
        int <- object$info$formula$vc[[pid]]$intercept
        if (inherits(object$info$model, "olmm")) {
            if (int == "ce") {
                parm[[pid]] <- c(paste("Eta", 1:object$info$model$dims["nEta"],
                                       ":(Intercept)", sep = ""), parm[[pid]])
            } else if (int == "ge") {
                parm[[pid]] <- c("(Intercept)", parm)
            }
        } else {
            if (int != "none") parm[[pid]] <- c("(Intercept)", parm[[pid]])
        }
    }
    return(parm)
}


##'-------------------------------------------------------- #
##' Extracts the estimates a fitted 'tvcm' object and
##' creates a list used in various methods.
##' variances of a fitted 'tvcm' object.
##'
##' @param object an object of class 'partynode'.
##' @param what   the type of estimate to be extracted.
##'
##' @return A list of 'fe' (fixed effects), 'vc' (varying
##' coefficients) and 're' (random effects) coefficients. The
##' 'vc' slot consists of separate matrices for each varying
##' coefficient partition.
##'
##' @details Used in 'coef', 'extract'.
##'-------------------------------------------------------- #

tvcm_get_estimates <- function(object, what = c("coef", "sd", "var"), ...) {
  
  what <- match.arg(what)
  model <- extract(object, "model")

  rval <- list(fe = numeric(),
               vc = replicate(length(object$info$node), matrix(,0,0)),
               re = numeric())
  names(rval$vc) <-  LETTERS[seq_along(object$info$node)]
  
  ## extract coefficients
  estimates <- switch(what,
                      coef = coef(model),
                      sd = diag(vcov(model)),
                      var = diag(vcov(model)))
  
  ids <- lapply(object$info$node, nodeids, terminal = TRUE)

  formList <- object$info$formula
    
  termsC <- tvcm_get_terms(names(coef(model)), ids, object$info$control$parm)
  termsE <- tvcm_get_terms(names(estimates), ids, object$info$control$parm)
    
  ## restricted coefficients
  if (any(termsE$type == "fe"))
    rval$fe <- estimates[termsE$type == "fe"]

  ## random effects
  if (any(termsE$type == "re"))
    rval$re <- estimates[termsE$type == "re"]

  ## varying coefficients
  if (any(termsE$type == "vc")) {

    for (pid in seq_along(object$info$node)) {

      nnodes <- length(ids[[pid]]) # number of nodes
      
      ## extract the terms corresponding to the partition
      vcTermsC <- unique(termsC$terms[termsC$partition == LETTERS[pid]])
      vcTermsE <- unique(termsE$terms[termsE$partition == LETTERS[pid]])

      ## build a matrix to store the coefficients
      rval$vc[[pid]] <- matrix(, nnodes, length(vcTermsC))
      rownames(rval$vc[[pid]]) <- ids[[pid]]
      colnames(rval$vc[[pid]]) <- vcTermsC

      ## fill the matrix
      for (i in seq_along(vcTermsC)) {
        subs <- termsE$terms == vcTermsC[i] & termsE$partition == LETTERS[pid]
        rval$vc[[pid]][termsE$node[subs], i] <- estimates[subs]
      }

      ## add colnames for varying intercepts
      subs <- which(colnames(rval$vc[[pid]]) %in% "")
      if (length(subs) > 0L) colnames(rval$vc[[pid]])[subs] <- "(Intercept)"
      if (ncol(rval$vc[[pid]]) > 0 & inherits(object$info$family, "family.olmm")) {
        tmp <- strsplit(colnames(rval$vc[[pid]]), ":")
        subs <- sapply(tmp, length) == 1L &
          sapply(tmp, function(x) substr(x[1L], 1L, 3L) == "Eta")
        colnames(rval$vc[[pid]])[subs] <-
          paste(colnames(rval$vc[[pid]])[subs], "(Intercept)", sep = ":")
      }

      ## fill the last row if necessary
      subs <- is.na(rval$vc[[pid]][nnodes, ])
      if (any(subs)) {
        
        ## compute the coefficients of the omitted node
        con <- model$contrasts[[paste("Node", LETTERS[pid],
                                      sep = "")]][nnodes, ]
        for (i in which(subs)) {
          rval$vc[[pid]][nnodes, i] <-
            switch(what,
                   coef = sum(con * rval$vc[[pid]][-nnodes, i], na.rm = TRUE),
                   sd = sum(con^2 * rval$vc[[pid]][-nnodes, i], na.rm = TRUE),
                   var = sum(con^2 * rval$vc[[pid]][-nnodes, i], na.rm = TRUE))
        }
      }
    }
    
    if (what == "sd") {
      getSqrt <- function(x) {
        if (is.list(x)) {
          return(lapply(x, getSqrt))
        } else if (!is.null(x)) {
          if (length(x) == 0) return(NA) else return(sqrt(x))
        } else {
          return(x)
        }
      }
      rval <- lapply(rval, getSqrt)
      
    }
}
  return(rval)
}


##'-------------------------------------------------------- #
##' Create short labes for 'vc' terms
##'
##' @param object a \code{\link{tvcm}} object
##' 
##' @return A character vector with a label for each 'vc'
##'    term.
##'
##' @details Used in 'tvcm_print' and 'plot.tvcm'.
##'-------------------------------------------------------- #

tvcm_print_vclabs <- function(object) {
  
  formula <- object$info$formula
  if (length(formula$vc) == 0) return(NULL)
  
  ## conditioning variables
  cond <- lapply(formula$vc, function(x) {
    rval <- all.vars(x$cond)
    if (length(rval) > 2L) rval <- c(rval[1L], "...")
    return(rval)
  })
  
  ## 'by' terms
  vcLabs <- terms(formula$original, specials = "vc")
  if (length(attr(vcLabs, "specials")$vc) == 0L) return(NULL)
  vcLabs <- rownames(attr(vcLabs, "factors"))[attr(vcLabs, "specials")$vc]
  vcLabs <- paste("getBy", vcLabs, sep = "_")
  getBy_vc <- function(..., by, intercept, nuisance) {
    mc <- match.call()
    rval <- deparse(mc$by, 500L)
    if (rval == "NULL") return("") else return(rval)
  }
  by <- sapply(vcLabs, function(x) eval(parse(text = x)))

  ## collapse the short labels
  rval <- rep.int("vc(", length(formula$vc))
  for (pid in seq_along(rval)) {
    if (length(cond) > 0L) 
      rval[pid] <- paste(rval[pid], paste(cond[[pid]], collapse = ", "), sep = "")
    if (by[pid] != "")
      rval[pid] <- paste(rval[pid], ", by = ", by[pid], sep = "")
  }
  rval <- paste(rval, ")", sep = "")
  return(rval)
}


##'-------------------------------------------------------- #
##' The main function to prune an object of class 'partynode'
##' according to some criteria.
##'
##' @param nodes   an object of class 'partynode'.
##' @param alpha   criteria according to which 'nodes' is to be
##'    pruned.
##' @param maxstep the maximal number of step (splits) to be
##'    accomplished
##' 
##' @return A list of class 'partynode'.
##'
##' @details Used in 'tvcm.prune'.
##'-------------------------------------------------------- #

tvcm_prune_node <- function(object, alpha = NULL, maxstep = NULL, terminal = NULL) {

  stopifnot(class(object)[1] %in% c("tvcm", "party", "partynode"))
  
  if ("partynode" %in% class(object)) {
    rval <- list(object)
  } else {
    rval <- object$info$node    
  }

  if (all(c(is.null(alpha), is.null(maxstep), is.null(terminal))))
    return(object$info$node)

  control <- extract(object, "control")

  if (!is.null(alpha) && depth(rval[[pid]]) > 0L) {
      splitpath <- object$info$splitpath
      p.value <- extract(object, "p.value")
      ms <- c(0, which(!is.na(p.value) & p.value <= alpha))
      if (length(ms) > 0L)
          ms <- max(ms[c(1, diff(ms)) == 1])
      maxstep <- min(maxstep, ms)
  }
  
  for (pid in seq_along(rval)) {

    ## prune the tree
    if (!is.null(maxstep))
      rval[[pid]] <- tvcm_prune_maxstep(rval[[pid]], maxstep)

    if (!is.null(terminal[[pid]]))
      rval[[pid]] <- tvcm_prune_terminal(rval[[pid]], terminal[[pid]])
    
    ## delete empty nodes and adjust id labeling
    tmp <- as.list(rval[[pid]])
    tmp <- tmp[sapply(tmp, function(x) !is.null(x))] # delete nodes
    ids_old <- unlist(lapply(tmp, function(x) x$id))
    ids_new <- 1L:length(tmp)
    for (i in ids_new) {
        tmp[[i]]$info$id$last <- tmp[[i]]$id
        tmp[[i]]$id <- ids_new[ids_old == tmp[[i]]$id]
        for (j in 1:length(tmp[[i]]$kids))
            tmp[[i]]$kids[j] <- ids_new[ids_old == tmp[[i]]$kids[j]]
    }
    rval[[pid]] <- as.partynode(tmp)
  }
  
  return(rval)
}


##'-------------------------------------------------------- #
##' Recursive function to delete nodes according to
##' the step in the algorithm the nodes were created.
##'
##' @param node    an object of class 'partynode'.
##' @param maxstep an integer scalar. The maximum step allowed.
##' 
##' @return A list of class 'partynode'.
##'
##' @details Used in 'tvcm_prune_node'.
##'-------------------------------------------------------- #

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


##'-------------------------------------------------------- #
##' Recursive function to delete nodes according node labels
##'
##' @param node     an object of class 'partynode'.
##' @param terminal an integer vector. The nodes considered
##'    as terminal nodes
##' 
##' @return A list of class 'partynode'.
##'
##' @details Used in 'tvcm_prune_node'.
##'-------------------------------------------------------- #

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
    }  }
  return(node)
}


##'-------------------------------------------------------- #
##' Creates a 'splitpath.tvcm' object for tracing the
##' fitting process.
##'
##' @param splitpath the splitpath obtained by the fitting process.
##' @param nodes     an object of class 'partynode'.
##' @param partData  a 'data.frame' with the partitioning variables.
##' @param control   an object of class 'tvcm_control'.
##'
##' @return A list of class 'splitpath.tvcm'.
##'
##' @details Used in 'tvcm' and 'prune.tvcm'.
##'-------------------------------------------------------- #

tvcm_grow_splitpath <- function(splitpath, varid, nodes, partData, control) {
    
  if (all(lapply(nodes, width) < 2L)) return(splitpath)  

  ## get the terminal node labels for each varying coefficient partition
  ids <- lapply(nodes, nodeids)
  
  ## assign each terminal node the step in which it was created
  steps <- lapply(nodes, function(node) {
    ids <- nodeids(node)
    rval <- nodeapply(node, ids, function(node) node$split$info$step)
    rval <- sapply(rval, function(x) if (is.null(x)) Inf else x)
    return(rval)
  })
  
  ## upate each step
  for (step in seq_along(splitpath)) {
    
    ## the selected partition
    partid <- splitpath[[step]]$partid
    
    ## get the nodes at disposition at the current step
    kidids <- vector("list", length(nodes))
    for (pid in seq_along(nodes)) {
      parentids <- ids[[pid]][steps[[pid]] < step]
      if (length(parentids) == 0L) {
        kidids[[pid]] <- 1L
      } else {
        kids <- nodeapply(nodes[[pid]], parentids, function(node) node$kids)
        for (i in seq_along(kids))
          kidids[[pid]] <-
            c(kidids[[pid]], unlist(lapply(kids[[i]], function(x) x$id)))
        kidids[[pid]] <- setdiff(sort(kidids[[pid]]), parentids)
      }
    }          
    
    ## re-label the descriptions of splits
    if (!is.null(splitpath[[step]]$varid)) {
      splitpath[[step]]$node <- kidids[[partid]][splitpath[[step]]$nodeid]
      splitpath[[step]]$var <- colnames(partData)[splitpath[[step]]$varid]
      splitpath[[step]]$cutpoint <- character_split(nodeapply(nodes[[partid]], splitpath[[step]]$node, function(node) node$split)[[1]], data = partData)$levels
    }
    
    if (control$sctest) {

      ## change row names of p-values tables
      for (pid in seq_along(nodes))
        if (!is.null(splitpath[[step]]$sctest[[pid]]))
          rownames(splitpath[[step]]$sctest[[pid]]) <-
            paste("Node", LETTERS[pid], kidids[[pid]], sep = "")
      
    }

    
    ## change the names of the loss grid elements !!!
    if (!is.null(splitpath[[step]]$lossgrid)) {
      names(splitpath[[step]]$lossgrid) <- LETTERS[seq_along(nodes)]
      for (pid in seq_along(splitpath[[step]]$lossgrid)) {
        names(splitpath[[step]]$lossgrid[[pid]]) <-
          paste("Node", kidids[[pid]], sep = "")
        for (nid in seq_along(splitpath[[step]]$lossgrid[[pid]])) { 
          names(splitpath[[step]]$lossgrid[[pid]][[nid]]) <-
            colnames(partData)[varid[[pid]]]
        }
      }
    }
  }
    
  class(splitpath) <- "splitpath.tvcm"
  return(splitpath)
}
