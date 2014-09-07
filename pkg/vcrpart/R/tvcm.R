##' -------------------------------------------------------- #
##' Author:      Reto Buergin
##' E-Mail:      reto.buergin@unige.ch, rbuergin@gmx.ch
##' Date:        2014-09-06
##'
##' Description:
##' The 'tvcm' function
##'
##' tvcolmm      convenience function for 'tvcm'
##' tvcglm       convenience function for 'tvcm'
##' tvcm         the main fitting function
##' tvcm_control control function for 'tvcm'
##'
##' all functions are documented as *.Rd files
##'
##' Last modifications:
##' 2014-09-06: - incorporate automatic cross-validation and pruning
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
  mcall <- call(name = fit,
               formula = quote(ff$full),
               family = quote(family),
               data = quote(mf),
               weights = weights,
               start = quote(start))
  mce <- match.call(expand.dots = TRUE)
  dotargs <- setdiff(names(mce), names(mc))
  dotargs <- intersect(dotargs, names(formals(fit)))
  dotargs <- setdiff(dotargs, names(mcall))
  for (arg in dotargs) mcall[[arg]] <- mce[[arg]]
  environment(mcall) <- environment()
  
  ## call root model
  model <- tvcm_grow_fit(mcall, doFit = FALSE)
  
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

  ## grow the tree

  ## replicates the required structure of 'tvcm' objects
  object <- structure(list(data = partData,
                           info = list(
                             call = mc,
                             mcall = mcall,                             
                             formula = formList,
                             direct = direct,
                             fit = fit,
                             family = family,
                             control = control,                             
                             model = model,
                             dotargs = dotargs)),
                      class = "tvcm")

  if (control$cv) {
    if (control$verbose) cat("\n* starting partitioning and cross validation ...\n")
    cvCall <- call(name = "cvloss",
                   object = quote(object),
                   folds = quote(control$folds),
                   type = "loss",
                   original = TRUE, 
                   verbose = FALSE,
                   papply = quote(control$papply))
    papplyArgs <- intersect(names(formals(control$papply)), names(control))
    papplyArgs <- setdiff(papplyArgs, names(args))
    for (arg in papplyArgs) cvCall[[arg]] <- control[[arg]]

    ## calls cvloss
    tree <- eval(cvCall)
    if (control$verbose)
      cat("\nestimated dfsplit =", format(tree$info$cv$dfsplit.hat, digits = 3), "\n")
    
  } else {
    
    ## calls directly 'tvcm_grow'
    if (control$verbose) cat("\n* starting partitioning ...\n")
    tree <- tvcm_grow(object)
    if (control$verbose) cat("OK")
  }
   
  ## pruning
  if (control$prune && inherits(tree, "tvcm")) {
    if (control$verbose) cat("\n* pruning ... ")
    tree <- prune(tree, dfsplit = tree$info$cv$dfsplit.hat, papply = control$papply)
    if (control$verbose) cat("OK")
  }
  
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
                         dfpar = 2.0, dfsplit = 0.0,
                         maxoverstep = ifelse(sctest, Inf, 0),
                         sctest = FALSE, alpha = 0.05, bonferroni = TRUE,
                         trim = 0.1, estfun = list(), ninpute = 5L,
                         maxfacsplit = 5L, maxordsplit = 10, maxnumsplit = 10,
                         cv = !sctest, folds = folds_control("kfold", 5),
                         prune = cv, keeploss = FALSE, papply = mclapply,
                         verbose = FALSE, ...) {
  mc <- match.call()
  
  ## check available arguments
  stopifnot(is.function(lossfun))
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
  stopifnot(is.logical(cv) && length(cv) == 1L)
  stopifnot(inherits(folds, "folds"))
  stopifnot(is.logical(prune) && length(prune) == 1L)
  if (!cv & prune) stop("'prune = TRUE' requires 'cv = TRUE'")
  stopifnot((is.logical(keeploss) | is.numeric(keeploss)) &&
            length(keeploss) == 1L)
  keeploss <- as.numeric(keeploss)
  if (is.numeric(verbose)) verbose <- as.logical(verbose)
  stopifnot(is.logical(verbose) && length(verbose) == 1L)
  
  ## check hidden arguments
  ptry <- ifelse(is.null(list(...)$ptry), Inf,  list(...)$ptry)
  stopifnot(is.numeric(ptry) && length(ptry) == 1L && ptry > 0)
  ntry <- if (is.null(list(...)$ntry)) Inf else list(...)$ntry
  stopifnot(is.numeric(ntry) && all(ntry > 0))
  vtry <- if (is.null(list(...)$vtry)) Inf else list(...)$vtry
  stopifnot(is.numeric(vtry) && all(vtry > 0))

  ## check and set 'papply'
  stopifnot(is.character(papply) | is.function(papply))
  if (is.function(papply)) {
    if ("papply" %in% names(mc)) {
      papply <- deparse(mc$papply)
    } else {
      papply <- deparse(formals(tvcm_control)$papply)
    }
  }
  
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
                  cv = cv,
                  folds = folds,
                  prune = prune,
                  keeploss = keeploss,
                  papply = papply,
                  verbose = verbose,
                  parm = NULL, intercept = NULL, 
                  functional.factor = "LMuo",
                  functional.ordered = "LMuo",
                  functional.numeric = "supLM")),
          class = "tvcm_control"))
}
