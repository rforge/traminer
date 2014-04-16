## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2013-12-05
##
## Description:
## Plot functions for 'tvcolmm' objects.
##
## References:
## party:           http://CRAN.R-project.org/package=party
## partykit:        http://CRAN.R-project.org/package=partykit
##
## Overview:
## plot.tvcolmm:   generic plot for tree.tvcolmm objects
## tvcolmm_terms_plot:  partial coefficient plots
## panel_get_main:
## panel_coef:
## panel_empty:
## --------------------------------------------------------- #

plot.tvcolmm <- function(x, type = c("default", "coef", 
                              "simple", "terms"),
                         main = NULL, drop_terminal = TRUE,
                         tnex = 1, newpage = TRUE,
                         pop = TRUE, gp = gpar(),
                         terminal = FALSE, 
                         legend = TRUE, legend_gp = gpar(),
                         terms = NULL, variable = NULL, ...) {
  
  type <- match.arg(type)
  
  if (type == "terms") {
    
    fun <- tvcolmm_terms_plot
    args <- list(object = x, terms = terms,
                 variable = variable)
    args <- append(args, list(...))
    do.call("fun", args)
  } else {
    
    ## tree plots
    
    ## terminal panel
    tp_fun <-
      switch(type,
             "default" = if (depth(x) > 3L) panel_empty else panel_coef,
             "coef" = panel_coef,
             "simple" = panel_empty)
    tp_args <-
      list(...)[names(list(...)) %in% names(formals(tp_fun))[-1]]
    if (type == "coef" | type == "default" && depth(x) <= 3L)
      tp_args$terms <- terms
      
    ## inner panel
    ip_fun <-
      switch(type,
             "default" = node_inner,
             "coef" = node_inner,
             "simple" = node_inner)
    ip_args <-
      list(...)[names(list(...)) %in% names(formals(ip_fun))[-1]]
    
    ## other arguments
    dot_args <-
      list(...)[names(list(...)) %in% names(formals(plot.party))[-1]]
    
    ## prepare call
    args <- list(x = x, main = main,
                 terminal_panel = tp_fun, tp_args = tp_args,
                 inner_panel = ip_fun, ip_args = ip_args,
                 drop_terminal = drop_terminal, tnex = tnex,
                 newpage = newpage,
                 pop = pop, gp = gp)
    args <- append(args, dot_args)
    
    ## call
    do.call("plot.party", args)
  }
}

tvcolmm_terms_plot <- function(object, terms = NULL,
                               variable = NULL, ask = NULL,
                               prob = NULL, neval = 50L, add = FALSE, ...) {

  data <- object$data
  allTerms <- names(coef(extract(object, "model")))
  allVars <- extract(object, "selected")
  allVars <- intersect(colnames(data), allVars)
  
  ## plot all varying coefficients ...
  partTerms <- grep("Part[1-9]?[0-9]?[0-9]", allTerms, value = TRUE)
  restTerms <- setdiff(allTerms, partTerms)
  partTerms <- unique(sub("Part[1-9]?[0-9]?[0-9]", "", partTerms))
  partTerms <- sub(":", "", partTerms)
  int <- grep("Eta[1-9]?[0-9]?[0-9]", partTerms)
  partTerms[int] <- paste(partTerms[int], ":Part", sep = "")
  partTerms[partTerms == ""] <- "Part"
  allTerms <- c(restTerms, partTerms)
  
  if (is.null(terms)) terms <- partTerms

  ## check terms
  if (any(!terms %in% allTerms)) stop(paste("terms '", paste(terms[!terms %in% allTerms], collapse = "', '"), "' are invalid. Use some of '", paste(allTerms, collapse = "', '"), "'.", sep = ""))
  
  ## ... and all partitioning variables that appear in the tree
  if (is.null(variable)) variable <- allVars
    
  ## set defaults
  if (is.null(ask)) ask <-
    ifelse(length(terms) * length(variable) == 1L, FALSE, TRUE)
  
  ## dotlist
  dotList <- list(...)

  ## get value spaces
  space <- tvcolmm_value_space(data, neval)

  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  for (i in 1:length(variable)) {

    z <- space[[variable[i]]]

    ## draw a random subsample to save time
    if (is.null(prob)) {
      if (length(z) * nrow(data) > 10000) {
        probi <- 10000 / (length(z) * nrow(data))
      } else {
        probi <- 1.0
      }
    } else {
      probi <- prob
    }
    tmp <- if (probi < 1) { data[tvcolmm_folds(object, K = 1L, type = "subsampling", prob = probi) == 1L, , drop = FALSE] } else { data }

    newdata <- NULL
    for (j in 1:length(z)) newdata <- rbind(newdata, tmp)
    newdata[, variable[i]] <- rep(z, each = nrow(tmp))
    beta <- predict(object, newdata = newdata, type = "terms")
    
    ## for each coefficients ...
    for (j in 1:length(terms)) {

      ## .. compute the expectation of beta_i at a value z_j
      ## over the distribution  of z_-j
      betaj <- tapply(beta[, terms[j]], newdata[, variable[i]], mean)

      ## define arguments for plot or points call
      ylab <-
        paste("coef(", terms[j], "|", variable[i], ")", sep = "")

      if (is.numeric(z)) {
        args <- appendDefArgs(dotList, list(type = "s"))
        args <- appendDefArgs(list(x = z, y = betaj), args)
      } else {
        args <- dotList
        args <- appendDefArgs(list(x = as.integer(z), y = betaj,
                                            type = "b"), args)
        if (!add) args <- appendDefArgs(args, list(axes = FALSE))
      }

      if (add) {
        fun <- points
      } else {
        fun <- plot.default
      }
      
      args <- appendDefArgs(args, list(xlab = variable[i], ylab = ylab))

      ## call plot or points
      do.call(fun, args)
      if (!add && is.factor(z) && !"axes" %in% names(list(...)) &&
          !args$axes) {
        box()
        axis(1, 1:nlevels(z), levels(z))
        axis(2)
      }
    }
  }
}

panel_get_main <- function(obj, node, id, dims) {
  rval <- ""
  if (id) rval <-
    paste(rval, paste(names(obj)[id_node(node)], sep = ""))
  if (id && any(c("n", "N") %in% dims))
    rval <- paste(rval, ": ", sep = "")
  if (all(c("n", "N") %in% dims)) {
    rval <- paste(rval,
                  "nobs = ", node$info$dims["n"],
                  "/ N = ", node$info$dims["N"], sep = "")
  } else if ("n" %in% dims) {
    rval <- paste(rval, "nobs = ", node$info$dims["n"])
  } else if ("N" %in% dims) {
    rval <- paste(rval, "N = ", node$info$dims["N"])
  }
  return(rval)
}

panel_coef <- function(obj, terms = NULL, id = TRUE,
                       dims = c("n", "N"),
                       margins = c(3, 2, 0, 0), gp = gpar(),
                       mean = FALSE, conf.int = FALSE, plot_gp = list(),
                       mean_gp = list(), conf.int_gp = list(), ...) {

  gp_def <- gpar(cex = 0.75)
  gp <- appendDefArgs(gp, gp_def)
  class(gp) <- "gpar"

  ## extract coefficients
  coef <- extract(obj, "coef")$varying
  if (is.null(terms)) terms <- colnames(coef)
  if (!is.list(terms)) terms <- list(terms)
  coef <- coef[, unlist(terms), drop = FALSE]
  coefList <- lapply(terms, function(trms) coef[, trms, drop = FALSE])
  
  if (conf.int) {
    sd <- extract(obj, "sd")$varying
    sd <- sd[, unlist(terms), drop = FALSE]
    sdList <- lapply(terms, function(trms) sd[, trms, drop = FALSE])
  }
  
  ## plot parameters
  if (length(plot_gp) == 0L) {
    plot_gp <- vector("list", length(terms))
  } else if (any(!unlist(lapply(plot_gp, is.list)))) {
    plot_gp <- lapply(1:length(terms), function(i) plot_gp)
  }
  
  plot_gp <-
    lapply(1:length(coefList),
           function(i) {
             if (conf.int) {
               ylim <- range(c(c(coefList[[i]] - 2 * sdList[[i]]),
                               c(coefList[[i]] + 2 * sdList[[i]])))
             } else {
               ylim <- range(coefList[[i]])
             }
             ylim <- ylim + c(-1, 1) * 0.1 * range(ylim)
             rval <- list(xlim = c(0.75, ncol(coefList[[i]]) + 0.25),
                          pch = 4L, ylim = ylim,
                          ylab = "coef",  type = "p",
                          label = abbreviate(colnames(coefList[[i]])),
                          height = 1, width = 0.6,
                          gp = gpar(cex = 1L))
             rval$ylim <- rval$ylim + 0.1 * c(-1, 1) * diff(rval$ylim)
             rval$ylim <- range(c(0, rval$ylim))
             if (length(plot_gp) > 0L)
               rval <- appendDefArgs(plot_gp[[i]], rval)
             return(rval)
           })

  if (conf.int) {

    if (length(conf.int_gp) == 0L) {
      conf.int_gp <- vector("list", length(terms))
    } else if (any(!unlist(lapply(conf.int_gp, is.list)))) {
      conf.int_gp <- lapply(1:length(terms), function(i) conf.int_gp)
    }    
    conf.int_gp_def <- list(angle = 90,
                            length = unit(2, "strwidth", "-"),
                            ends = "both",
                            type = "open")
    conf.int_gp <- lapply(1:length(terms),
                          function(i) appendDefArgs(conf.int_gp[[i]],
                                                    conf.int_gp_def))
  }
  
  ## population mean
  if (mean) {

    if (length(mean_gp) == 0L) {
      mean_gp <- vector("list", length(terms))
    } else if (any(!unlist(lapply(mean_gp, is.list)))) {
      mean_gp <- lapply(1:length(terms), function(i) mean_gp)
    }    
    mean_gp_def <- list(gp = gpar(col = "grey50", lty = 1, lwd = 0.75), pch = 4L)
    mean_gp <- lapply(1:length(terms),
                      function(i) appendDefArgs(mean_gp[[i]], mean_gp_def))

    w <- tapply(weights(obj), predict(obj, type = "node"), sum) / sum(weights(obj))
    meanCoef <- colSums(coef * matrix(w, nrow(coef), ncol(coef)))
    meanCoef <- lapply(terms, function(trms) meanCoef[trms, drop = FALSE])
    if (conf.int) {
      meanSd <- sqrt(colSums(sd^2 * matrix(w, nrow(coef), ncol(coef))^2))
      meanSd <- lapply(terms, function(trms) sqrt(meanSd[trms, drop = FALSE]))
    }
  }

  qN <- qnorm(0.975)
  
  rval <- function(node) {
    
    strUnit <- unit(2, "strheight", "A")
    pushViewport(viewport(layout = grid.layout(2, 1, heights = unit.c(strUnit, unit(1, "npc") - strUnit)))) # 1
    grid.rect(gp = gpar(fill = "white", col = 0))
    
    pushViewport(viewport(layout.pos.row = 1L)) # 2
    grid.text(panel_get_main(obj, node, id, dims))
    upViewport()
    
    pushViewport(viewport(layout.pos.row = 2L)) # 2
    pushViewport(viewport(layout = grid.layout(length(coefList), 1))) # 3
    
    for (i in 1:length(coefList)) {
      
      pushViewport(viewport(layout.pos.row = i)) # 4
      
      pushViewport(viewport(height = unit(plot_gp[[i]]$height, "npc"),
                            width = unit(plot_gp[[i]]$width, "npc"))) # 5
      pushViewport(plotViewport(margins = margins,
                                xscale = plot_gp[[i]]$xlim, yscale = plot_gp[[i]]$ylim,
                                default.units = "native")) # 6
      
      grid.segments(unit(1:ncol(coefList[[i]]), "native"),
                    unit(rep(plot_gp[[i]]$ylim[1], ncol(coefList[[i]])), "native"),
                    unit(1:ncol(coefList[[i]]), "native"),
                    unit(rep(plot_gp[[i]]$ylim[2], ncol(coefList[[i]])), "native"),
                    gp = gpar(col = "lightgrey"))
      
      grid.segments(unit(plot_gp[[i]]$xlim[1], "native"), unit(0, "native"),
                    unit(plot_gp[[i]]$xlim[2], "native"), unit(0, "native"),
                    gp = gpar(col = "black"))
      
      grid.segments(unit(plot_gp[[i]]$xlim[1], "native"), unit(0, "native"),
                    unit(plot_gp[[i]]$xlim[2], "native"), unit(0, "native"),
                    gp = gpar(col = "black"))

      if (conf.int) {
        
        grid.segments(unit(1:ncol(coefList[[i]]), "native"),
                      unit(coefList[[i]][as.character(id_node(node)),] -
                           qN * sdList[[i]][as.character(id_node(node)),], "native"),
                      unit(1:ncol(coefList[[i]]), "native"),
                      unit(coefList[[i]][as.character(id_node(node)),] +
                           qN * sdList[[i]][as.character(id_node(node)),], "native"),
                      arrow = arrow(angle = conf.int_gp[[i]]$angle,
                        length = conf.int_gp[[i]]$length, 
                        ends = conf.int_gp[[i]]$ends,
                        type = conf.int_gp[[i]]$type),
                      gp = plot_gp[[i]]$gp)
        
        if (FALSE) {
          
          grid.segments(unit(1:ncol(coefList[[i]]), "native"),
                        unit(meanCoef[[i]] - qN * meanSd[[i]], "native"),
                        unit(1:ncol(coefList[[i]]), "native"),
                        unit(meanCoef[[i]] + qN * meanSd[[i]], "native"),
                        arrow = arrow(angle = conf.int_gp[[1]]$angle,
                        length = conf.int_gp[[i]]$length, 
                        ends = conf.int_gp[[i]]$ends,
                        type = conf.int_gp[[i]]$type),
                        gp = mean_gp[[i]]$gp)

        }
        
      }
      
      if (plot_gp[[i]]$type %in% c("p", "b")) {

        if (mean) {
          grid.points(unit(1:ncol(coefList[[i]]), "native"),
                      unit(meanCoef[[i]], "native"),
                      pch = mean_gp[[i]]$pch, gp = mean_gp[[i]]$gp)
        }
        
        grid.points(unit(1:ncol(coefList[[i]]), "native"),
                    unit(coefList[[i]][as.character(id_node(node)),],
                         "native"),
                    pch = plot_gp[[i]]$pch, gp = plot_gp[[i]]$gp)
      }

      
      if (plot_gp[[i]]$type %in% c("l", "b")) {

        if (mean) {
          grid.lines(unit(1:ncol(coefList[[i]]), "native"),
                     unit(meanCoef[[i]], "native"),
                     gp = mean_gp[[i]]$gp)
        }
        
        grid.lines(unit(1:ncol(coefList[[i]]), "native"),
                   unit(coefList[[i]][as.character(id_node(node)),],
                        "native"),
                   gp = plot_gp[[i]]$gp)
      }
      
      if (id_node(node) == min(nodeids(obj, terminal = TRUE))) {
        grid.yaxis(gp = gpar(lineheight = 0.5))
        grid.text(plot_gp[[i]]$ylab, unit(-2, "lines"), unit(0.5, "npc"),
                  rot = 90)
      }
      
      grid.xaxis(at = 1:ncol(coefList[[i]]), label = plot_gp[[i]]$label,
                 gp = gpar(lineheight = 0.5))
      
      grid.rect()
      upViewport(3L)

    }
    
    ## close viewports
    upViewport(3L)
    return(rval)
  }
}
class(panel_coef) <- "grapcon_generator"


panel_empty <- function(obj, id = TRUE, dims = c("n", "N"), ...) {
    
  rval <- function(node) {
    
    nid <- id_node(node)
    pushViewport(viewport(height = unit(3, "strheight", "A")))
    grid.rect(gp = gpar(fill = "white", col = 0))
    grid.text(panel_get_main(obj, node, id, dims))
    upViewport()
  }
  return(rval)
}
class(panel_empty) <- "grapcon_generator"
