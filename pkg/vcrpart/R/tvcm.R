## --------------------------------------------------------- #
## Author:      Reto Buergin
## E-Mail:      reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:        2014-05-03
##
## Description:
## The tvcm function
##
## tvcm:         the main fitting function
## tvcm_control: control function for 'tvcm'
## --------------------------------------------------------- #

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
  
  ## check arguments
  if (control$verbose) cat("* checking arguments ... ")
  mc <- match.call(expand.dots = FALSE)
  stopifnot(inherits(formula, "formula"))
  stopifnot(inherits(control, "tvcm_control"))
  env <- environment(eval.parent(mc$formula))

  ## set fitting function
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
    
  ## set family
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
  if (any(substr(all.vars(formula), 1, 4) == "Node"))
    stop("'Node' is a reserved variable name and cannot be used as",
         "variable name nor as prefix of a variable name (sorry).")
  formList <- vcrpart_formula(formula, family, env)
  ff <- tvcm_formula(formList, family, env)
  
  ## extract model frames
  if (control$verbose) cat("OK\n* extracting model frames ... ")
  m <- match(c("data", "subset", "weights", "na.action"), names(mc), 0L)
  mf <- mc[c(1L, m)]
  mf$formula <- formList$all
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")  
  mf <- eval.parent(mf)
  
  ## create a call
  start <- list(...)$start
  weights <- model.weights(mf)
  if (is.null(weights)) weights <- rep(1.0, nrow(mf))  
  call <- call(name = fit,
               formula = quote(ff$root),
               family = quote(family),
               data = quote(mf),
               weights = weights,
               start = quote(start))
  mce <- match.call(expand.dots = TRUE)
  dotargs <- setdiff(names(mce), names(mc))
  for (arg in dotargs) call[[arg]] <- mce[[arg]]
  environment(call) <- environment()
  
  ## call root model
  if (inherits(family, "family.olmm")) call$doFit <- FALSE
  model <- eval(call)
  if (inherits(family, "family.olmm")) call$doFit <- TRUE  

  ## get coefficients to be used in the root node test
  control <- tvcm_update_control(control, model, formList)
  
  ## define model data
  mf <- mf[rownames(model.frame(model)),, drop = FALSE]

  ## define partitioning data
  partForm <- if (!is.null(formList$vc)) formList$vc$cond else formula(~ 1)
  partvar <- model.frame(partForm, mf)
  attr(partvar, "terms") <- attr(mf, "terms")
  attr(partvar, "na.action") <- attr(mf, "na.action")
  
  if (control$verbose) cat("OK\n* starting partitioning ...\n")
  
  ## set the root node
  nodes <- partynode(id = 1L, info = list(dims = nobs(model), depth = 0L))
  mf$Node <- factor(rep(1L, nrow(mf)))
  splitpath <- list()
  
  run <- 1L
  step <- 0L
  
  while (run > 0L) {

      step <- step + 1L

      test <- NULL
      split <- NULL
      
      if (control$verbose) cat("\n* starting step", step, "...")
      
      ## get current partitions
      where <- fitted_node(nodes, partvar)
      mf$Node <- factor(where)
      
      varid <- 1L:ncol(partvar)
      nodeid <- 1L:width(nodes)
   
      ## set start values if required
      if (control$fast > 0L)
        start <- tvcm_get_start(model, levels(mf$Node), control$parm)
      
      ## --------------------------------------------------- #
      ## Step 1: fit the current model
      ## --------------------------------------------------- #

      model <- tvcm_fit_model(call, control)
      if (width(nodes) == control$maxwidth | step > control$maxstep |
          length(unique(sapply(splitpath, function(x) x$varid))) == control$nselect |
          control$maxdepth == 0L | length(control$parm) == 0L) { 
          run <- 0L
      }

      call$formula <- quote(ff$tree)

      if (run > 0L && control$method == "mob") { # you can set this to false

        ## --------------------------------------------------- #
        ## Step 2: variable selection via coefficient constancy tests
        ## --------------------------------------------------- #
        
        test <- tvcm_fit_sctest(model, nodes, partvar, control, call)
        
        ## return error if test failed
        if (inherits(test, "try-error"))
          stop("coefficient constancy tests failed. Abort")
        
        run <- 1L * (min(c(1+.Machine$double.eps, test$p.value),
                         na.rm = TRUE) <= control$alpha)
          
        if (run > 0L) {
          
          pval <-  apply(test$p.value, 2,
                         function(x) suppressWarnings(min(x, na.rm = TRUE)))
          varid <- if (length(pval) > 0) which.min(pval) else NULL
          nodeid <- if (length(pval) > 0) which.min(test$p.value[, varid]) else NULL
        } 
      }
      
      
      if (run > 0L) {
        
        ## ------------------------------------------------- #
        ## Step 3: search a cutpoint
        ## ------------------------------------------------- #
        
        split <- tvcm_fit_risk(varid, nodeid, model, nodes, partvar, 
                                  control, call, step)

        risk <- unlist(lapply(split$riskgrid,
                              lapply, function(x) c(Inf, x[, ncol(x)])))
        if (!any(na.omit(risk) < Inf)) run <- 0L
      }
      
      
      if (run > 0L)
        nodes <- tvcm_fit_splitnode(nodes, split, partvar, step)

      if (run >= 0L)
          splitpath[[step]] <-
            list(step = step,                          
                 varid = split$varid,
                 nodeid = split$nodeid,
                 cutid = split$cutid,
                 sctest = test,
                 logLik = logLik(model),
                 risk = control$riskfun(model),
                 riskgrid = split$riskgrid)
      
      ## print the actions
      if (control$verbose) {
        if (run > 0L) {
          ids <- nodeids(nodes)
          cat("\n\nSplitting variable:", names(partvar)[split$varid])
          cat("\nNode:", ids[split$nodeid])
          cat("\nCutpoint:\n\n")
          print(split$riskgrid[[split$varid]][[split$nodeid]][split$cutid,,drop=0])
        } else {
          cat("\nNo admissible split found. Return object.\n")
        }
      }
  }
  
  ## if 'intercept == "po"', refit the model with appropriate contrasts
  if (nlevels(mf$Node) > 1L && formList$vc$intercept == "ge") {
    con <- contr.sum(levels(mf$Node))
    tab <- tapply(weights, mf$Node, sum)
    con[nrow(con),] <- con[nrow(con),] * tab[-length(tab)] / tab[length(tab)]
    colnames(con) <- levels(mf$Node)[1:(nlevels(mf$Node) - 1)]
    contrasts <- eval(call$contrasts)
    contrasts$Node <- con
    call$contrasts <- contrasts
    model <- tvcm_fit_model(call, control)
  }
  
  if (control$verbose) cat("\n* building object ...")

  ## inscribe original node names for later pruning
  nodes <- as.list(nodes)
  for (i in 1:length(nodes)) nodes[[i]]$info$original$id <- nodes[[i]]$id
  nodes <- as.partynode(nodes)
  
  ## prepare the title
  title <- c("Tree-based varying-coefficients model")
  
  ## modify splitpath    
  splitpath <- tvcm_update_splitpath(splitpath, nodes, partvar, control)
  
  ## the output object
  tree <- party(nodes, data = partvar,
                fitted = data.frame(
                  "(fitted)" = fitted_node(nodes, data = partvar),
                  "(response)" = model.response(model.frame(model)),
                  "(weights)" = weights(model),
                  check.names = FALSE),
                terms = terms(formula, keep.prder = TRUE),
                info = list(
                  title = title,
                  call = mc,
                  formula = ff,
                  fit = fit,
                  family = family,
                  control = control,
                  model = model,
                  nstep = width(nodes) - 1L,
                  splitpath = splitpath,
                  dotargs = list(...)))   
  class(tree) <- c("tvcm", "party")

  if (control$verbose) {
    cat("OK\n\nFitted model:\n")
    print(tree)
  }
  
  if (control$verbose)
    cat("* computations finished, return object\n")
  
  return(tree)
}

tvcm_control <- function(method = c("mob", "partreg"),
                         alpha = 0.05, bonferroni = TRUE,
                         maxwidth = ifelse(method == "partreg", 10L, Inf),
                         minbucket = 50L, mtry = Inf,
                         minsplit = 2 * minbucket, trim = 0.1,
                         maxdepth = Inf, maxstep = maxwidth - 1L, 
                         nselect = Inf, estfun = list(),  
                         maxevalsplit = 20L, riskfun = deviance,
                         fast = 0L, verbose = FALSE, ...) {

  method <- match.arg(method)
    
  ## check available arguments
  stopifnot(alpha >= 0.0 & alpha <= 1.0)
  stopifnot(is.logical(bonferroni))
  stopifnot(minsplit >= 0)

  ## set the default parameters for 'gefp.estfun' calls
  estfun <- appendDefArgs(estfun, list(predecor = TRUE,
                                       nuisance = NULL,
                                       silent = TRUE))
  if (!is.null(estfun$level) && estfun$level != "observation")
    warning("'level' argument for 'estfun' is set to 'observation'")
  estfun$level <- "observation"

  ## create a list of parameters of class 'tvcm_control'
  return(structure(
           appendDefArgs(
             list(...),
             list(method = method,
                  alpha = alpha,
                  bonferroni = bonferroni,
                  maxwidth = maxwidth,
                  minbucket = minbucket,
                  mtry = mtry,
                  minsplit = minsplit,                 
                  trim = trim,
                  maxdepth = maxdepth,
                  maxstep = maxstep,
                  nselect = nselect,
                  estfun = estfun,
                  maxevalsplit = maxevalsplit,
                  riskfun = riskfun,
                  fast = fast, 
                  verbose = verbose,
                  parm = NULL, intercept = NULL,
                  functional.factor = "LMuo",
                  functional.ordered = "LMuo",
                  functional.numeric = "supLM")),
          class = "tvcm_control"))
}
