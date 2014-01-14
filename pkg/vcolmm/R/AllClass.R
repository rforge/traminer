## class for ordinal linear mixed models fitted with 'olmm'
setClass("olmm",
         representation(
           env = "environment",
           frame = "data.frame",     # model frame
           call = "call",            # the user call
           formula = "formula",      # model formula
           terms = "list",           # terms of X and W
           y = "ordered",            # response vector
           X = "matrix",             # fixef mm
           W = "matrix",             # ranef model matrix
           subject = "factor",       # grouping factor
           subjectName = "character",# variable name of 'subject'
           weights = "numeric",      # weights
           weights_sbj = "numeric",  # weights on group level
           offset = "numeric",       # offset
           xlevels = "list",         # levels of cat. covariates
           contrasts = "list",       # contrasts
           dims = "integer",         # dimensions
           fixef = "matrix",         # odds-variable effects
           ranefCholFac = "matrix",  # cholesky factor matrix
           coefficients = "numeric", # model parameters
           restricted = "logical",   # restricted coefficients
           eta = "matrix",           # linear predictor
           u = "matrix",             # stand. random effects
           logLik_obs = "numeric",   # observation-wise ll
           logLik_sbj = "numeric",   # subject-wise ll
           logLik = "numeric",       # marginal log Likelihood
           score_obs = "matrix",     # observation-wise scores
           score_sbj = "matrix",     # subject-wise scores.
           score = "numeric",        # (total) sum of scores
           info = "matrix",          # expected fisher info
           ghx = "matrix",           # zeros used for GHQ
           ghw = "matrix",           # weights used for GHQ
           ranefElMat = "matrix",    # eliminaton matrix
           optim = "list",           # call for optimization
           output = "list"           # output of optimisation
           ))

## class for summaries on 'olmm' objects
setClass("summary.olmm",
         representation(
           methTitle = "character",
           formula = "character",
           data = "character",
           subset = "character",
           AICtab = "data.frame",
           FEmatEtaInv = "matrix",
           FEmatEtaVar = "matrix",
           REmat = "matrix",
           dims = "integer"))
