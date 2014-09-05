##' -------------------------------------------------------- #
##' Author:      Reto Buergin
##' E-Mail:      reto.buergin@unige.ch, rbuergin@gmx.ch
##' Date:        2014-09-02
##'
##' Description:
##' The 'tvcm' function
##'
##' tvcolmm      convenience function for 'tvcm'
##' tvcglm       convenience function for 'tvcm'
##' tvcm         the main fitting function
##' tvcm_control control function for 'tvcm'
##'
##' Last modifications:
##' 2014-09-04: - assign only those arguments of '...' to 'fit'
##'               that appear in 'formals(fit)'
##' 2014-08-02: - the 'formula' slot is now a list of formulas as
##'               produced by 'vcrpart_formula'. The modification
##'               was due to acceleration techniques ('vcrpart_formula'
##'               is usually slow!)
##' 2014-08-29: - implement adjustment of loss reduction by number
##'               of predictor of coefficient-group
##' 2014-07-31: - set 'sctest = FALSE' as the default
##'             - return an error if multiple trees and 'sctest = TRUE'
##'             - check if global intercept is removed if
##'               intercepts are tested with coef. const. tests
##'             - add new arguments 'dfsplit' and 'maxoverstep'
##'               to 'tvcm_control'
##'             - add new stopping criteria based on 'dfsplit'
##'               and 'maxoverstep'
##' 2014-06-26: incorporate new function 'tvcm_grow_setsplits'
##' 2014-06-16: allow coefficient-wise trees
##' -------------------------------------------------------- #

tvcolmm <- function(formula, data, family = cumulative(),
                    weights, subset, na.action,
                    control = tvcm_control(), ...) {
    mc <- match.call()
    mc[[1L]] <- as.name("tvcm")
    if (!"family" %in% names(mc) &
        (length(mc) < 4L |
         length(mc) >= 4L && !inherits(eval.parent(mc[[4L]]), "family.olmm")))
        mc$family <- formals(tvcolmm)$family        
    mc$fit <- "olmm"
    return(eval.parent(mc))
}

tvcglm <- function(formula, data, family,
                   weights, subset, na.action,
                   control = tvcm_control(), ...) { 
    mc <- match.call()
    mc[[1L]] <- as.name("tvcm")       
    mc$fit <- "glm"
    return(eval.parent(mc))
}

tvcm <- function(formula, data, fit, family, 
                 weights, subset, na.action,
                 control = tvcm_control(), ...) {

  ## get specified arguments
  mc <- match.call(expand.dots = FALSE)

  ## check and set arguments  
  if (control$verbose) cat("* checking arguments ... ")
  stopifnot(inherits(formula, "formula"))
  stopifnot(inherits(control, "tvcm_control"))
  
  ## check and set 'fit'
  if (missing(fit)) {
    if (missing(family)) stop("no 'family'.")
    fit <- switch(class(family),
                  family.olmm = "olmm",
                  family = "glm",
                  stop("no 'fit' function"))
  } else {
    if (is.function(fit)) fit <- deparse(mc$fit)
  }
  if (!fit %in% c("glm", "olmm")) stop("'fit' not recognized.")
    
  ## check and set 'family'
  if (missing(family)) stop("no 'family'.")
  if (is.character(family)) {
      family <- get(family, mode = "function", envir = parent.frame())
  } else if (is.function(family)) {
      family <- family()
  }
  if (!class(family) %in% c("family", "family.olmm")) stop("'family' not recognized")
  if (fit != "olmm") control$estfun <- NULL
  
  ## set formulas
  if (control$verbose) cat("OK\n* setting formulas ... ")
  if (any(grepl("Right", all.vars(formula)) | grepl("Left", all.vars(formula)) |
          grepl("Node", all.vars(formula))))
  if (any(substr(all.vars(formula), 1, 4) == "Node"))
    stop("'Node', 'Left' and 'Right' are reserved labeles and cannot be used as",
         "substrings of variable names.")
  env <- environment(eval.parent(mc$formula))
  formList <- vcrpart_formula(formula, family, env)
  nPart <- length(formList$vc)  

  direct <- any(sapply(formList$vc, function(x) x$direct))
  if (length(direct) == 0L) direct <- FALSE
  control$direct <- direct

  vcRoot <- rep(TRUE, nPart)
  ff <- tvcm_formula(formList, vcRoot, family, env)
  
  ## extract model frames
  if (control$verbose) cat("OK\n* extracting model frames ... ")
  m <- match(c("data", "subset", "weights", "na.action"), names(mc), 0L)
  mf <- mc[c(1L, m)]
  mf$formula <- formList$all
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")  
  mf <- eval.parent(mf)

  y <- as.data.frame(model.response(mf))
  if (ncol(y) > 1L) mf <- cbind(y, mf[, -1, drop = FALSE])
  
  ## create a call
  if (control$verbose) cat("OK\n* setting arguments ... ")
  start <- list(...)$start
  weights <- model.weights(mf)
  if (is.null(weights)) weights <- rep(1.0, nrow(mf))  
  call <- call(name = fit,
               formula = quote(ff$full),
               family = quote(family),
               data = quote(mf),
               weights = weights,
               start = quote(start))
  mce <- match.call(expand.dots = TRUE)
  dotargs <- setdiff(names(mce), names(mc))
  dotargs <- intersect(dotargs, names(formals(fit)))
  for (arg in dotargs) call[[arg]] <- mce[[arg]]
  environment(call) <- environment()
  
  ## call root model
  model <- tvcm_fit_model(call, doFit = FALSE)
  
  ## check if there are categorical variables among the predictors
  etaVars <- unlist(lapply(formList$vc, function(x) {
    lapply(x$eta, function(x) all.vars(x))
  }))
  etaVars <- intersect(etaVars, colnames(model.frame(model))) 
  if (any(sapply(model.frame(model)[, etaVars, drop = FALSE], is.factor)))
    stop("variables in 'by' of 'vc' terms must be numeric. ",
         "Use 'model.matrix' to convert the categorical variables to ",
         "numeric predictors.")

  ## set whether coefficient constancy tests are used    
  if (control$sctest) {
    if (nPart > 1L)
      stop("coefficient constancy tests can be used only ",
           "if a single 'vc' term is specified.")
    if (!is.null(formList$vc) && (formList$vc[[1L]]$direct & formList$fe$intercept != "none"))
      stop("if 'sctest = TRUE', searching for intercept is only possible if the",
           "global intercept is removed. Use something like 'formula = y ~ -1 + ...'")
  }
  
  ## set 'parm' for the root node
  control <- tvcm_grow_setcontrol(control, model, formList, vcRoot, FALSE)

  ## set imputation in 'control'
  if (!inherits(model, "olmm") | inherits(model, "olmm") &&
      length(unique(table(model$subject))) == 1L)
    control$ninpute <- 1L
  
  ## specify which coefficients are considered as 'nuisance' parameters
  if (control$verbose && control$sctest)
    if (length(control$estfun$nuisance) > 0L)
      cat("\n\tnuisance parameters: ",
          paste(paste("'", control$estfun$nuisance, "'", sep = ""),
                collapse = ", ", sep = ""))

  if (control$verbose && !control$sctest)
    if (length(unlist(control$nuisance)) > 0L)
      cat("\n\tnuisance terms:",
          paste(paste("'", unlist(control$nuisance), "'", sep = ""),
                collapse = ", ", sep = ""))
  
  ## define model data
  mf <- mf[rownames(model.frame(model)),, drop = FALSE]
  
  ## partitioning variables
  partVars <- lapply(formList$vc, function(x) attr(terms(x$cond), "term.labels"))
  
  ## define partitioning data
  if (!is.null(formList$vc)) {
    partForm <- formula(paste("~", paste(unique(unlist(partVars)), collapse = "+")))
  } else {  
    partForm <- formula(~ 1)
  }
  partData <- model.frame(partForm, mf)
  if (any(sapply(partData, function(x) !(is.factor(x) | is.numeric(x)))))
    stop("partitioning variables must be either 'numeric' or (ordered) 'factor'.")
  
  attr(partData, "terms") <- attr(mf, "terms")
  attr(partData, "na.action") <- attr(mf, "na.action")
  varid <- lapply(partVars, function(x) {
    as.integer(sapply(x, function(x) which(colnames(partData) == x))) })
  
  if (control$verbose) cat("\n* starting partitioning ...\n")
  
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
      where[[pid]] <- vcrpart_contr.sum(where[[pid]], weights(model))
      mf[, paste("Node", LETTERS[pid], sep = "")] <- where[[pid]]
    }
    
    nodeid <- lapply(nodes, function(x) 1:width(x))
    
    if (control$verbose) cat("\n* starting step", step, "...")
    
    ## --------------------------------------------------- #
    ## Step 1: fit the current model
    ## --------------------------------------------------- #

    vcRoot <- sapply(nodeid, length) == 1L
    ff <- tvcm_formula(formList, vcRoot, family, env)
    model <- try(tvcm_fit_model(call))
    
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
                                 control, call, formList, step), silent = TRUE)
      
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
          noverstep <- if (loss$loss -
                           control$dfpar * loss$df -
                           control$dfsplit * 1 < 0)
            noverstep + 1L else 0L
          if (noverstep > control$maxoverstep) {
            run <- 0L
            stopinfo <- "'maxoverstep' reached"
            nstep <- nstep - 1L
          }
        }
        
      }
    }
    
    ## incorporate the split into 'nodes'
    if (run > 0L)
      nodes <- tvcm_grow_splitnode(nodes, where, loss, partData,
                                   step, weights(model))

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
                         row.names = paste("step", step + c(-1, 0))))
        
      } else {
        cat("\n\nStopping the algorithm.\nMessage:", as.character(stopinfo), "\n")
        if (inherits("try-error", stopinfo)) warning(as.character(stopinfo))
        
      }
    }
  }
  
  if (control$verbose) cat("\n* building object ... ")
  
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
    tree <- party(nodes[[1L]],
                  data = partData,
                  fitted = data.frame(
                    "(response)" = model.response(model.frame(model)),
                    "(weights)" = weights(model),
                    check.names = FALSE),
                  terms = terms(formula, keep.order = TRUE),
                  info = list(
                    title = title,
                    call = mc,
                    formula = formList,
                    direct = direct,
                    fit = fit,
                    family = family,
                    control = control,
                    info = stopinfo,
                    model = model,
                    node = nodes,
                    nstep = nstep,
                    splitpath = splitpath,
                    pruned = FALSE,
                    dotargs = list(...)))
    class(tree) <- c("tvcm", "party")
  }
  if (control$verbose) cat("OK")
  
  if (control$verbose) {
    cat("\n\nFitted model:\n")
    print(tree)
  }
  
  if (control$verbose)
    cat("* computations finished, return object\n")
  
  return(tree)
}

tvcm_control <- function(lossfun = neglogLik2, 
                         maxstep = Inf, maxwidth = Inf,
                         minsize = 30, maxdepth = Inf,
                         dfpar = 2.0, dfsplit = 0.0, maxoverstep = 0,
                         sctest = FALSE, alpha = 0.05, bonferroni = TRUE,
                         trim = 0.1, estfun = list(), ninpute = 5L,
                         maxfacsplit = 5L, maxordsplit = 10, maxnumsplit = 10,
                         keeploss = FALSE, verbose = FALSE, ...) {
  
  ## check available arguments
  stopifnot(is.function(lossfun))
  stopifnot((is.logical(keeploss) | is.numeric(keeploss)) &&
            length(keeploss) == 1L)
  keeploss <- as.numeric(keeploss)
  stopifnot(is.numeric(maxstep) && length(maxstep) == 1L && maxstep >= 0L)
  stopifnot(is.numeric(maxwidth) && all(maxwidth > 0L))
  stopifnot(is.null(minsize) | (is.numeric(minsize) && all(minsize > 0)))
  stopifnot(is.numeric(maxdepth) &&  all(maxdepth >= 0))
  stopifnot(is.numeric(dfpar) && length(dfpar) == 1L)
  stopifnot(is.numeric(dfsplit) && length(dfsplit) == 1L)
  stopifnot(is.numeric(maxoverstep) && length(maxoverstep) == 1L)
  stopifnot(is.logical(sctest) && length(sctest) == 1L)
  if (length(alpha) != 1L)
  stopifnot(is.numeric(alpha) && length(alpha) == 1L && alpha >= 0.0 && alpha <= 1.0)
  stopifnot(is.logical(bonferroni) && length(bonferroni) == 1L)
  stopifnot(is.numeric(trim) && length(trim) == 1L && trim >= 0.0 & trim < 0.5)
  stopifnot(is.list(estfun))
  stopifnot(is.numeric(ninpute) && length(ninpute) == 1L)
  stopifnot(is.numeric(maxfacsplit) && length(maxfacsplit) == 1L && maxfacsplit > 1L)
  stopifnot(is.numeric(maxordsplit) && length(maxordsplit) == 1L && maxordsplit > 1L)
  stopifnot(is.numeric(maxnumsplit) && length(maxnumsplit) == 1L && maxnumsplit > 1L)
  if (is.numeric(verbose)) verbose <- as.logical(verbose)
  stopifnot(is.logical(verbose) && length(verbose) == 1L)
  
  ## check hidden arguments
  ptry <- ifelse(is.null(list(...)$ptry), Inf,  list(...)$ptry)
  stopifnot(is.numeric(ptry) && length(ptry) == 1L && ptry > 0)
  ntry <- if (is.null(list(...)$ntry)) Inf else list(...)$ntry
  stopifnot(is.numeric(ntry) && all(ntry > 0))
  vtry <- if (is.null(list(...)$vtry)) Inf else list(...)$vtry
  stopifnot(is.numeric(vtry) && all(vtry > 0))
  
  ## set the default parameters for 'gefp.estfun' calls
  estfun <- appendDefArgs(estfun, list(predecor = TRUE,
                                       nuisance = NULL,
                                       silent = FALSE))
  if (!is.null(estfun$level) && estfun$level != "observation")
    warning("'level' argument for 'estfun' is set to 'observation'")
  estfun$level <- "observation"

  ## ensure backward compability
  if ("maxevalsplit" %in% names(list(...))) maxnumsplit <- list(...)$maxevalsplit
  if ("minbucket" %in% names(list(...))) minsize <- list(...)$minbucket
  
  ## create a list of parameters of class 'tvcm_control'
  return(structure(
           appendDefArgs(
             list(...),
             list(lossfun = lossfun,
                  keeploss = keeploss,
                  maxstep = maxstep,
                  maxwidth = maxwidth,
                  minsize = minsize,
                  maxdepth = maxdepth,
                  dfpar = dfpar,
                  dfsplit = dfsplit,
                  maxoverstep = maxoverstep,
                  ptry = ptry,
                  ntry = ntry,
                  vtry = vtry,
                  trim = trim,
                  sctest = sctest,
                  alpha = alpha,
                  bonferroni = bonferroni,
                  estfun = estfun,
                  ninpute = ninpute,
                  maxfacsplit = maxfacsplit,
                  maxordsplit = maxordsplit,
                  maxnumsplit = maxnumsplit,
                  verbose = verbose,
                  parm = NULL, intercept = NULL, 
                  functional.factor = "LMuo",
                  functional.ordered = "LMuo",
                  functional.numeric = "supLM")),
          class = "tvcm_control"))
}
