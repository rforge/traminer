## From Joint Sequence Analysis, (R. Piccarreta, SMR, 2017)
##  measures of association between dimensions
##  portions of code inspired from the assoc.domains function of seqhandbook

dissdomassoc <- function(domdiss, jointdiss = NULL, assoc = c("pearson","R2"),
      names=NULL, weights=NULL) {

## domdiss: list of dissimilarities matrices or objects (one per channel)
## jointdiss: dissimilarity matrix or object between sequences of combined states
## assoc: requested measure of association

  if (!is.list(domdiss) || length(domdiss) < 2)
    stop("domdiss should be a list of at least two distance matrices or objects")

  if (!all(is.numeric(unlist(domdiss, use.names=FALSE))))
    stop("non numeric values in distance matrices")

  #assoclist <- c("pearson","spearman","kendall","R2","cronbach","cron.subsets","all")
  ## kendall takes too much time for large number of diss values
  assoclist <- c("pearson","spearman","R2","cronbach","cron.subsets","all")
  if (!all(assoc %in% assoclist))
    stop("bad assoc values, allowed are ", paste(assoclist, collapse=","))
  if ("all" %in% assoc) assoc <- assoclist[c(1,2,3,5)]
  if ("R2" %in% assoc & !any(c("pearson","spearman") %in% assoc) ) {
    stop("R2 can only be used in combination with 'pearson' or 'spearman'")
  }

  ## setting weights
  if (inherits(domdiss[[1]],'dist'))
    ncases <- attr(domdiss[[1]],'Size')
  else
    ncases <- nrow(domdiss[[1]])
  if (is.null(weights)) {
    weights <- rep(1, ncases)
    weighted=FALSE
    }
  else {
    if (length(weights) != ncases)
      stop("length of weights not equal to number of cases!")
    weighted=TRUE
  }
  ww <- as.numeric(as.dist(weights %*% t(weights)))
  sww <- sqrt(ww)

  ndom <- length(domdiss)
  ndomv <- ndom - 1
  ## transforming into vector of distances
  distlist <- lapply(domdiss, function(x) as.numeric(as.dist(x)))
  dissmat <- matrix(unlist(distlist, use.names=FALSE), ncol=ndom, byrow=FALSE)
  if (is.null(names)) names <- names(domdiss)
  colnames(dissmat) <- names

  if (!is.null(jointdiss)){
    if (!all(is.numeric(jointdiss)))
      stop("when not NULL, jointdiss must be a distance matrix or object")
    jointdiss <- as.numeric(as.dist(jointdiss))
    dissmat <- cbind(dissmat,jointdiss)
    ndomv <- c(rep(ndom-1,ndom),ndom)
  }
  dissmat <- apply(dissmat,2,scale)

  if ("spearman" %in% assoc) { ## we replace columns with weighted ranked
    rankmat <- apply(dissmat,2,rank)
    if (weighted) {
      dissmat.spear <- apply(dissmat,2,weighted.rank)
      ## weighted.rank returns NA for min and max ranks
      ## we replace these NAs with the non-weighted ranks
      dissmat.spear[is.na(dissmat.spear)] <- rankmat[is.na(dissmat.spear)]
    } else {
      dissmat.spear <- rankmat
    }
    #rm(rankmat)
    dissmat.spear <- apply(dissmat.spear,2,scale)
    ##dissmat.spearw <- sww * dissmat.spear
  }
  ##dissmatw <- sww * dissmat

  ## defining function
  rsquare.corr <- function(correlation, jointdiss, ndom, ndomv) {
    corr.tmp <- correlation
    if (!is.null(jointdiss)) {
      corr.tmp[1:ndom,ndom+1] <- 0
    }
    return((rowSums(corr.tmp^2)-1)/ndomv)
  }

  res <- list()
  if ("pearson" %in% assoc){ ## cor.wt from the psych package
    ##correlation <- cor(dissmatw, method='pearson')
    correlation <- cor.wt(dissmat, w=ww)$r
    res[["Pearson"]] <- correlation
    if ("R2" %in% assoc) res[["Pearson.Rsquare"]] <- rsquare.corr(correlation, jointdiss, ndom, ndomv)
  }
  if ("spearman" %in% assoc){
    ##correlation <- cor(dissmat.spearw)
    correlation <- cor.wt(dissmat.spear, w=ww)$r
    res[["Spearman"]] <- correlation
    if ("R2" %in% assoc) res[["Spearman.Rsquare"]] <- rsquare.corr(correlation, jointdiss, ndom, ndomv)
  }
  #if ("kendall" %in% assoc){  ## kendall takes too much time
  #  correlation <- cor(dissmat, method='kendall')
  #  res[["Kendall"]] <- correlation
  #  if ("R2" %in% assoc) res[["Kendall.Rsquare"]] <- rsquare.corr(correlation, jointdiss, ndom, ndomv)
  #}
  if (any(c("cronbach","cron.subsets") %in% assoc)){
    ## Cronbach alpha for all domains
    #sdissmat <- scale(dissmat)
    sigmatot <- var(rowSums(dissmat[,1:ndom]))
    chron <- (ndom/(ndom-1))*(1-ndom/sigmatot)
    names(chron) <- paste0('(',paste0(names,collapse=','),')')
    res[["Cronbach"]] <- chron
  }

  ## Cronbach alpha for each combination of domains
  if ("cron.subsets" %in% assoc) {
    cron.subset <- numeric()
    if (ndom>2) {
      for(p in (ndom-1):2) {
        set.start <- length(cron.subset) + 1
        sets <- combn(1:ndom, p, simplify=FALSE)
        set.end <- length(cron.subset) + length(sets)
        for(i in 1:length(sets)) {
          sigmatot <- var(rowSums(dissmat[,sets[[i]],drop=FALSE]))
          cronbach <- (p/(p-1))*(1-p/sigmatot)
          cron.subset <- c(cron.subset,cronbach)
          #names(cron.subset)[length(cron.subset)] <- paste0('(',paste0(names[sets[[i]]],collapse=','),')')
        }
        names(cron.subset)[set.start:set.end] <- lapply(sets, function(x) paste0('(',paste0(names[x],collapse=','),')'))
      }
      res[["Cronbach.subsets"]] <- cron.subset
    }
    else
      message("Two or less domains, no subset possible for Cronbach alpha")
  }

  #res <- list(correlation, Rsquare, res3)

  #names(res) <- c(corrname,'Rsquared',"Cronbach alpha")

  return(res)
}
