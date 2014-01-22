## --------------------------------------------------------- #
## Author:      Reto Buergin
## E-Mail:      reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:        2014-01-15
##
## Description:
## The tvcolmm function
##
## tvcolmm:         the main fitting function for ordinal linear
##                  mixed model trees with varying fixed effects
## tvcolmm_control: control function for 'tvcolmm'
##
## Todo:
## -
## --------------------------------------------------------- #

tvcolmm <- function(formula, data, control = tvcolmm_control(),
                    subset, na.action, offset, 
                    vi = c("none", "po", "npo"),
                    linear = NULL, ...) {

  ## check arguments
  if (control$verbose) cat("* checking arguments ... ")
  mc <- match.call(expand.dots = FALSE)
  stopifnot(inherits(formula, "formula"))
  stopifnot(inherits(control, "tvcolmm_control"))
  if (any(all.vars(formula) == "Part")) stop("'Part' is a reserved keyword and cannot be used for variable names")
  
  ## set restricted parameters
  control$restricted <- linear 
  
  ## set varying intercept
  vi <- match.arg(vi)
  control$intercept <- vi
  
  ## set formulas
  formula <- as.Formula(formula)
  if (length(formula)[2] < 2)
    stop("no splitting variables")

  ## set formulas
  if (control$verbose) cat("OK\n* setting formulas ... ")
  fullForm <- as.Formula(formula)
  formList <- olmm_formula(formula(fullForm, rhs = - length(fullForm)[2]),
                           env = parent.frame(n = 1))
  fullForm <- update(fullForm, formList$full)
  environment(fullForm) <- environment(eval.parent(mc$formula))
  partForm <- formula(as.Formula(fullForm),
                      lhs = -1, rhs = length(fullForm)[2])
  modelForm <- formula(as.Formula(formula),
                       lhs = 1, rhs = 1:(length(fullForm)[2] - 1))
  
  ## extract model frames
  if (control$verbose) cat("OK\n* extracting model frames ... ")
  m <- match(c("data", "subset", "weights", "na.action", "offset"),
             names(mc), 0)
  frame <- mc[c(1L, m)]
  frame$formula <- fullForm
  frame$drop.unused.levels <- TRUE
  frame[[1L]] <- as.name("model.frame")
  if (control$maxsurrogate > 0 | control$probsurrogate)
    frame$na.action <- na.pass  
  frame <- eval.parent(frame)

  ## compress category labels
  frame <- tvcolmm_modify_catpreds(formula = modelForm, data = frame)
    
  ## list with elements to fit the model
  args <- append(list(formula = modelForm, data = frame), list(...))
  
  ## call root model
  args$doFit <- FALSE
  model <- try(do.call("olmm", args))
  args$doFit <- TRUE
  args$weights <- weights(model)
      
  ## modify 'control' and arguments for fitting the model
  control <- tvcolmm_modify_control(model, control)
  args <- tvcolmm_modify_modargs(model, args)

  ## extract formulas
  ff <- list(original = formula,
             root = args$formula,
             tree = tvcolmm_formula(model, control))
  
  ## define model data
  frame <- frame[rownames(model.frame(model)), , drop = FALSE]

  ## define partitioning data
  partvar <- model.frame(partForm, frame)

  if (!is.null(control$type.vars) && # check
      (length(control$type.vars) != ncol(partvar) |
       any(!control$type.vars %in% c("tm-var", "tm-inv"))))
    stop("'type.vars' is missspecified")
  
  if (is.null(control$type.vars)) { # retrieve types if necessary
    getType <- function(i) {
      nval <- rowSums(table(model@subject, partvar[, i]) > 0)
      return(c("tm-inv", "tm-var")[1 + 1 * (max(nval) > 1)])
    }
    control$type.vars <- sapply(1:ncol(partvar), getType)
  }
  
  if (control$verbose) cat("\n* starting partitioning ...\n")
  
  ## set the root node
  nodes <- partynode(id = 1L,
                     info = list(dims =
                          c(N = nlevels(model@subject),
                            n = length(model@subject)),
                          depth = 0L))
  
  args$data$Part <- factor(rep(1L, nrow(args$data)))

  run <- TRUE
  step <- 0L
  
  while (run) {

      step <- step + 1L
      if (control$verbose) cat("\n* starting step", step, "...")
      
      ## get current partitions
      where <- fitted_node(nodes, partvar)
      args$data$Part <- factor(where)
      
      varid <- 1:ncol(partvar)
      partid <- 1:width(nodes)
   
      ## set start values if required
      if (control$fast > 0L)
        args$start <- tvcolmm_get_start(model, args)
      
      ## --------------------------------------------------- #
      ## Step 1: fit the current model
      ## --------------------------------------------------- #
      
      model <- tvcolmm_fit_model(ff, args, control, control$verbose)

      ## return error if fitting failed
      if (inherits(model, "try-error"))
        stop("model fitting failed.")

      ## --------------------------------------------------- #
      ## Step 2: variable selection via coefficient constancy tests
      ## --------------------------------------------------- #

      if (control$fluctest) {
        
        test <- tvcolmm_fit_fluctest(model, nodes, partvar, control,
                                     ff, args)
        
        ## return error if test failed
        if (inherits(test, "try-error"))
          stop("coefficient constancy tests failed.")
      
        run <- suppressWarnings(min(test$p.value, na.rm = TRUE)) <= control$alpha

        if (run) {

          pval <-  apply(test$p.value, 2, function(x) suppressWarnings(min(x, na.rm = TRUE)))
          varid <- which.min(pval)
          partid <- which.min(test$p.value[, varid])
        
          if (control$verbose) {

            cat("\nSplitting variable:", colnames(partvar)[varid])
            cat("\nPartition:", levels(args$data$Part)[partid])
            cat("\nSelection statistic (p-value) = ")
            cat(format(pval[varid], digits = 3))
          }
          
        }
      } else {

        test <- NULL
      }
      
      
      if (run) {
        
        ## ------------------------------------------------- #
        ## Step 3: search a cutpoint
        ## ------------------------------------------------- #
        
        newnodes <- tvcolmm_fit_splitnode(varid, partid,
                                          partvar, nodes,
                                          model, ff, args, test,
                                          control, step)
          
        if (is.partynode(newnodes)) {
          nodes <- newnodes
        
        } else {

          ## add fluctest in the info slot of each terminal node
          ids <- nodeids(nodes, terminal = TRUE)
          nodes <- as.list(nodes)
          for (i in 1:length(nodes)) {
            if (nodes[[i]]$id %in% ids) {
              nodes[[i]]$info$step <- step
              nodes[[i]]$info$fluctest <- test
            }
          }
          nodes <- as.partynode(nodes)

          ## set stop flag
          run <- FALSE
        }
      }

      ## check width of tree
      if (width(nodes) > control$maxwidth) {
        if (control$verbose)
          cat("\nMaximal width reached. Return object.\n")
        run <- FALSE
      }
    }

  ## if 'intercept == "po"', refit the model with appropriate contrasts
  if (nlevels(args$data$Part) > 0 && control$intercept == "po") {
    con <- contr.sum(levels(args$data$Part))
    tab <- tapply(args$weights, args$data$Part, sum)
    con[nrow(con),] <- con[nrow(con),] * tab[-length(tab)] / tab[length(tab)]
    colnames(con) <- levels(args$data$Part)[1:(nlevels(args$data$Part) - 1)]
    args$contrasts$Part <- con
    model <- tvcolmm_fit_model(ff, args, control, FALSE)
  }

  if (control$verbose) {
    cat("\nFitted model:\n")
    print(model)
  }

  if (control$verbose) cat("\n* building object ...")

  
  ## prepare the title
  title <-
    paste(c("Cumulative", "Baseline", "Adjacent")[model@dims["family"]])
  if (model@dims["family"] == 1)
    title <-
      paste(title, c("Logit", "Probit", "Cauchy")[model@dims["link"]])
  title <- paste(title, "Mixed Model with Varying-Effects")
  control$terms <- control$terms$original
  
  ## the output object
  tree <- party(nodes, data = partvar,
                fitted = data.frame(
                  "(fitted)" = fitted_node(nodes, data = partvar),
                  "(response)" = model@y,
                  "(weights)" = model@weights,
                  check.names = FALSE),
                terms = terms(as.Formula(formula), keep.order = TRUE),
                info = list(
                    title = title,
                    call = mc,
                    formula = ff,
                    vi = control$intercept,
                    linear = control$restricted,
                    control = control,
                    model = model,
                    fluctest = test,
                    dotargs = list(...)))

  class(tree) <- c("tvcolmm", "party")

  if (control$verbose)
    cat("OK\n* computations finished, return object\n")
  
  return(tree)
}

tvcolmm_control <- function(alpha = 0.05, bonferroni = TRUE,
                            minsplit = 50L, trim = 0.1,
                            lossfun = neglogLik, breakties = FALSE,
                            terms = NULL, verbose = FALSE,...) {
  
  ## check available arguments
  stopifnot(alpha >= 0 & alpha <= 1)
  stopifnot(is.logical(bonferroni))
  stopifnot(minsplit >= 0)
  
  rval <- appendDefArgs(
            list(...),
            list(alpha = alpha,
                 bonferroni = bonferroni,
                 minsplit = minsplit,
                 maxdepth = Inf,
                 maxwidth = Inf,
                 nselect = Inf,
                 maxevalsplit = 20,
                 lossfun = lossfun,
                 breakties = breakties,
                 terms = terms, restricted = NULL,
                 intercept = "none",
                 mtry = Inf,
                 maxsurrogate = 0L,
                 probsurrogate = FALSE,
                 condsurrogate = "observation",
                 functional.factor = "LMuo",
                 functional.ordered = "LMuo",
                 functional.numeric = "supLM",
                 type.vars = NULL,
                 fluctest = TRUE,
                 fast = 0L,
                 verbose = verbose))
  
  class(rval) <- "tvcolmm_control"
  return(rval)
}
