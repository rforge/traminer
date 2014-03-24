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
## panel_fluctest:
## --------------------------------------------------------- #

plot.tvcolmm <- function(x, type = c("default", "coef", 
                              "fluctest", "simple", "terms"),
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
    
    if (type == "default" & depth(x) < 1L)
      type <- "fluctest"
    
    ## tree plots
    
    ## terminal panel
    tp_fun <-
      switch(type,
             "default" = if (depth(x) > 3L) panel_empty else panel_coef,
             "coef" = panel_coef,
             "fluctest" = panel_fluctest,
             "simple" = panel_empty)
    tp_args <-
      list(...)[names(list(...)) %in% names(formals(tp_fun))[-1]]
    if (type == "coef" | type == "default" && depth(x) <= 3L)
      tp_args$terms <- terms
      
    ## inner panel
    ip_fun <-
      switch(type,
             "default" = if (depth(x) > 3L) node_inner else panel_fluctest,
             "coef" = node_inner,
             "fluctest" = if (terminal) node_inner else panel_fluctest,
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
                  "n = ", node$info$dims["n"],
                  "/ N = ", node$info$dims["N"], sep = "")
  } else if ("n" %in% dims) {
    rval <- paste(rval, "n = ", node$info$dims["n"])
  } else if ("N" %in% dims) {
    rval <- paste(rval, "N = ", node$info$dims["N"])
  }
  return(rval)
}

panel_coef <- function(obj, terms, id = TRUE, dims = c("n", "N"),
                       margins = c(2, 2, 0, 0), gp = gpar(), ...) {
    
    coef <- coef(obj, nodeids(obj))$varying[, terms, drop = FALSE]
    
    gp_def <- gpar(cex = 0.75)
    gp <- appendDefArgs(list(...), gp_def)
    class(gp) <- "gpar"
    
    args_def <- list(xlim = c(0.75, ncol(coef) + 0.25), pch = 4L,
                     ylim = range(c(coef)), ylab = "coef", 
                     label = abbreviate(colnames(coef)),
                     height = 1, width = 0.6)    
    args_def$ylim <- args_def$ylim + 0.1 * c(-1, 1) * diff(args_def$ylim)
    args_def$ylim <- range(c(0, args_def$ylim))
    args <- appendDefArgs(list(...), args_def)
    
    rval <- function(node) {
            
      strUnit <- unit(2, "strheight", "A")
      pushViewport(viewport(layout = grid.layout(2, 1, heights = unit.c(strUnit, unit(1, "npc") - strUnit))))
      grid.rect(gp = gpar(fill = "white", col = 0))
      
      pushViewport(viewport(layout.pos.row = 1L))
      grid.text(panel_get_main(obj, node, id, dims))
      upViewport()

      pushViewport(viewport(layout.pos.row = 2L))
      pushViewport(viewport(height = unit(args$height, "npc"),
                            width = unit(args$width, "npc")))
      pushViewport(plotViewport(margins = margins,
                                xscale = args$xlim, yscale = args$ylim,
                                default.units = "native"))

      grid.segments(unit(1:ncol(coef), "native"),
                    unit(rep(args$ylim[1], ncol(coef)), "native"),
                    unit(1:ncol(coef), "native"),
                    unit(rep(args$ylim[2], ncol(coef)), "native"),
                    gp = gpar(col = "lightgrey"))
      
      grid.segments(unit(args$xlim[1], "native"), unit(0, "native"),
                    unit(args$xlim[2], "native"), unit(0, "native"),
                    gp = gpar(col = "black"))
      
      grid.segments(unit(args$xlim[1], "native"), unit(0, "native"),
                    unit(args$xlim[2], "native"), unit(0, "native"),
                    gp = gpar(col = "black"))

      grid.points(unit(1:ncol(coef), "native"),
                  unit(coef[as.character(id_node(node)),],
                       "native"),
                  pch = args$pch, gp = gp)

      if (id_node(node) == min(nodeids(obj, terminal = TRUE))) {
          grid.yaxis(gp = gpar(lineheight = 0.5))
          grid.text(args$ylab, unit(-2, "lines"), unit(0.5, "npc"),
                    rot = 90)
      }

      grid.xaxis(at = 1:ncol(coef), label = args$label,
                 gp = gpar(lineheight = 0.5))

      grid.rect()
      
      ## close viewports
      upViewport(4L)
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

panel_fluctest <- function(obj, id = TRUE, dims = c("n", "N"),
                           margins = c(2, 2, 0, 2), log = FALSE,
                           all = NULL, gp = gpar(), ...) {
  
  if (is.null(all)) all <- depth(obj) < 3 & ncol(obj$data) < 5
  
  gp_def <- gpar(cex = 0.75, col = "black", fontface = 1L,
                 lwd = 1, fontsize = get.gpar()$fontsize)
  gp <- appendDefArgs(gp, gp_def)
  class(gp) <- "gpar"
  
  alpha <- extract(obj, "control")$alpha
  if (log) {
    if (alpha < 1e-3) alpha <- 1e-4
    alpha <- (log10(alpha) + 4) / 4
  }
  
  args_def <- list(height = 0.75, width = 0.75, pch = 4L)
  args_def$xlim <- if (all | log) c(0, 1) else c(0, 1.1 * alpha)
  
  args <- appendDefArgs(list(...), args_def)    
  FUN <- function(x) suppressWarnings(min(x, na.rm = TRUE))
  
  rval <- function(node) {
    
    pval <- order.by <- NULL
    testMessage <- ""
    show <- FALSE        
    if (is.terminal(node) && !is.null(obj$info$fluctest)) {
      pval <-
        obj$info$fluctest$p.value[paste("Part", node$id, sep = ""), ]
      order.by <- names(pval)
      show <- TRUE
    } else if (!is.null(split_node(node)$info$fluctest)) {
      pval <- apply(split_node(node)$info$fluctest$p.value, 2, FUN)
      order.by <- names(pval)
      show <- TRUE
    } else {
      testMessage <- "test not available"
    }
    
    if (show & log) {
      pval[pval < 1e-4] <- 1e-4
      pval <- (log10(pval) + 4) / 4
    }
    
    if (show & (all | (!all & any(pval <= alpha, na.rm = TRUE)))) {
      
      if (!all) {
        subs <- !is.na(pval) & pval <= alpha
        order.by <- order.by[subs]
        pval <- pval[subs]
      }
      best <- which.min(pval)
      if (!is.null(args$ylim)) ylim <- ylim else ylim <- c(0.5, length(pval) + 0.5)
      gp_node <- gp
      gp_node$cex <- rep(gp$cex, length(pval))
      gp_node$col <- rep(gp$col, length(pval))
      gp_node$lwd <- rep(gp$lwd, length(pval))
      gp_node$fontface <- rep(gp$fontface, length(pval))
      gp_node$fontsize <- rep(gp$fontsize, length(pval))
      pch <- rep(args$pch, length(pval))
      pch[is.na(pval)] <- 1L
      gp_node$lwd[best] <- 2 * max(gp$lwd)
      gp_node$fontface[best] <- 2L
      pval[is.na(pval)] <- args$xlim[2]
    } else {
      if (show) {
        show <- FALSE
        testMessage <- "min(p-value) > alpha"
      }
      ## pseudo values
      gp_node <- gp
      ylim <- c(0, 1)
    }
    
    pushViewport(viewport(
                   y = unit(0.5 + (args$height - 0.75) / 2, "npc"),
                   height = unit(args$height, "npc"),
                   width = unit(args$width, "npc")))
    strUnit <- unit(2, "strheight", "A")
    pushViewport(viewport(layout = grid.layout(2, 1, heights = unit.c(strUnit, unit(1, "npc") - strUnit))))
    grid.rect(gp = gpar(fill = "white", col = 0))
    
    pushViewport(viewport(layout.pos.row = 1))
    grid.text(panel_get_main(obj, node, id, dims))
    upViewport()
    
    pushViewport(viewport(layout.pos.row = 2))
    
    grid.rect(gp = gpar(fill = "white", col = 0))
    pushViewport(plotViewport(margins = margins, xscale = args$xlim,
                                  yscale = ylim, default.units = "native"))
    
    if (show) {
      
      grid.segments(unit(rep(args$xlim[1], length(pval)), "native"),
                    unit(1:length(pval), "native"),
                    unit(rep(args$xlim[2], length(pval)), "native"),
                    unit(1:length(pval), "native"),
                    gp = gpar(col = "lightgrey"))            
      grid.lines(unit(c(alpha, alpha), "native"),
                 unit(ylim, "native"),
                 gp = gpar(lty = 2))
      grid.points(unit(pval, "native"),
                  unit(1:length(pval), "native"),
                  pch = pch, gp = gp_node)
      grid.yaxis(at = 1:length(pval), label = order.by,
                 gp = gpar(lineheight = 0.5, 
                   fontsize = gp_node$fontsize,
                   fontface = gp_node$fontface))
      
      if (log) {
        at <- (c(-3, -2, -1, 0) + 4) / 4
        label <- c(0.001, 0.01, 0.1, 1)
        subs <- which(label > args$xlim[1] & label < args$xlim[2])
        grid.xaxis(at = at[subs],
                   label = label[subs],
                   gp = gpar(lineheight = 0.5))
      } else {                
        grid.xaxis(at = args$xlim, label = args$xlim,
                   gp = gpar(lineheight = 0.5))
      }
      grid.text(paste("p-value", if (log) "(log)", ""),
                y = unit(-1.5, "lines"))
      
      grid.rect()
      
    } else {          
      grid.rect(gp = gpar(fill = "white", col = 0))
      grid.text(testMessage,
                gp = gpar(fontface = 2,
                  fontsize = 1.1 * gp$fontsize))
    }
    
    upViewport(4L)
  }
  return(rval)
  }
class(panel_fluctest) <- "grapcon_generator"
