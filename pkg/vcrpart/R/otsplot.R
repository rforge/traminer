otsplot_control <- function(cex = 1, lwd = 1/4, col = NULL,
                            hide.col = grey(0.8),
                            lorder = c("background", "foreground") ,
                            lcourse = c("upwards", "downwards"),
                            grid.scale = 1/5, grid.lwd = 1/2,
                            grid.fill =  grey(0.95), grid.col = grey(0.6),          
                            layout = NULL, margins = c(5.1, 4.1, 4.1, 3.1),
                            strip.fontsize = 12, strip.fill =  grey(0.9),
                            pop = TRUE, newpage = TRUE, maxit = 500L) {

  ## match arguments with duplicated entries in defaults
  lorder <- match.arg(lorder)
  lcourse <- match.arg(lcourse)

  ## checks
  stopifnot(is.logical(pop))
  stopifnot(is.logical(newpage))
  stopifnot(is.null(layout) | (is.numeric(layout) && length(layout) == 2L))
  stopifnot(is.numeric(margins) && length(margins) == 4L)
  stopifnot(is.numeric(maxit) && maxit > 0L)

  ## set list for plot parameters
  plot_gp <- list(cex = cex, lwd = lwd, col = col,
                  hide = hide.col, noobs = "black",
                  lorder = lorder, lcourse = lcourse,
                  alpha = 1, sf_cex = 0.5, sf_cex_leaves = 1)
  
  ## set list for translation zone plot parameters
  grid_gp <- list(scale = grid.scale, lwd = grid.lwd,
                  fill = grid.fill, col = grid.col)
  
  ## set list for strips plot parameters
  strip_gp <- gpar(fontsize = strip.fontsize,
                   fill = strip.fill)
   
  ## return a 'otsplot_control' object
  return(structure(
           list(plot_gp = plot_gp,
                grid_gp = grid_gp,
                strip_gp = strip_gp,
                layout = layout,
                margins = margins,
                pop = pop,
                newpage = newpage,
                maxit = maxit),
           class = "otsplot_control"))
}

otsplot_filter <- function(method = c("minfreq", "cumfreq", "linear"), level = NULL) {
    method <- match.arg(method)
    if (is.null(level) && method %in% c("minfreq", "cumfreq"))
        stop("'otsplot_filter' requires an inpute for 'level'.")
    method <- switch(method,
                     linear = linear,
                     minfreq = minfreq,
                     cumfreq = cumfreq)
    return(structure(list(method = method, level = level), class = "otsplot_filter"))
}

otsplot.default <- function(x, y, subject, weights = NULL, groups,
                            control = otsplot_control(), filter = NULL, 
                            main, xlab, ylab, xlim, ylim, ...) {  

  mc <- match.call()
  
  ## check filter argument
  if (is.character(filter))
      filter <- otsplot_filter(method = filter[1L], level = list(...)$level)
  if (!is.null(filter))
      stopifnot(inherits(filter, "otsplot_filter"))
  stopifnot(inherits(control, "otsplot_control"))
  
  ## check control argument
  stopifnot(inherits(control, "otsplot_control"))
  cArgs <- list(...)
  cArgs <- cArgs[intersect(names(cArgs), names(formals(otsplot_control)))]
  cArgs <- do.call("olmm_control", cArgs)
  control <- appendDefArgs(cArgs, control)

  ## check compatibility between 'filter' and 'lorder' control argument
  if (!is.null(filter)) control$plot_gp$lorder <- "foreground"
  
  ## check and prepare raw data
  
  if (sum(c(length(x), length(y)) != length(subject)) != 0) {
    stop("input vectors 'x', 'y' and 'subject' have different lengths.")
  }
  
  ## check the x variable
  if (is.numeric(x)) {
    if (length(unique(x)) > 1) {
      xdiff <- diff(sort(unique(x)))
      if (sum((xdiff / min(xdiff) -
               round(xdiff / min(xdiff),0)) != 0) == 0) {
        x <- factor(x, levels = seq(min(x, na.rm = TRUE),
                         max(x, na.rm = TRUE), min(xdiff)))
      } else {
        warning("problems with distances between 'x' positions. ",
                "The 'x' positions will not be illustrated adequately.")
        x <- factor(x)
      }
    } else { x <- factor(x) }  # pseudo case
  }
  
  ## y must be an ordered factor
  if (!is.factor(y)) stop("'y' must be a 'factor'.") 
  
  ## convert subject variable to factor
  if (!is.factor(subject)) subject <- factor(subject)
  subject <- droplevels(subject)
  
  ## construct/ convert groups variable
  if (!missing(groups)) {
    if (!is.factor(groups)) groups <- factor(groups)
    groups <- droplevels(groups)
    if (sum(is.na(groups)) > 1) {
      groups <- factor(groups, levels = c(levels(groups), "not available"))
      groups[is.na(groups)] <- "not available"
    }
    if (length(groups) == length(subject)) 
      groups <- unique(data.frame(subject, groups))[, 2]    
    if (length(groups) != nlevels(subject)) 
      stop("cannot link 'groups' with 'subject'.")    
    
  } else groups <- factor(rep(1, nlevels(subject)))
  names(groups) <- levels(subject)
  
  ## construct/ convert weights variable
  if (!is.null(weights)) {
    if (!is.numeric(weights)) stop("'weights' are not numeric.")
    if (min(weights, na.rm = TRUE) < 0L) stop("negative 'weights'.")
    if (length(weights) == length(subject))
      weights <- unique(data.frame(subject, weights))[, 2L]
    if (length(weights) != nlevels(subject) | (sum(is.na(weights))) > 0L)
      stop(paste("cannot link 'weights' with 'subject'."))
  } else weights <- rep(1L, nlevels(subject))

  ## remove NAs and omit.levels in x, y or subject
  SUBS <- is.na(subject) | is.na(x) | is.na(y)
  if (sum(SUBS) > 0) {
    subject <- subject[!SUBS]
    x <- x[!SUBS]
    y <- y[!SUBS]
  }
  
  ## standardize weights
  wsubjectGroups.unscaled <- tapply(weights, groups, sum)
  TMP1 <- table(groups)[as.integer(groups)]
  TMP2 <- tapply(weights, groups, sum)[as.integer(groups)]
  weights <- weights / TMP2 * TMP1
  names(weights) <- levels(subject)
  
  ## subject variable
  subjectUsed <- intersect(levels(subject), unique(subject))
  nSubjectTot <- nlevels(subject)
  subject <- factor(subject, levels = subjectUsed)
  weightsUse <- weights[names(weights) %in% subjectUsed]
  nSubjectUse <- nlevels(subject)
  subject <- as.integer(subject)
  
  ## x variable
  xLevs <- levels(x)
  nXLevs <- nlevels(x)
  x <- as.integer(x)
  
  ## y variable
  yLevs <- levels(y)
  nYLevs <- nlevels(y)
  y <- as.integer(y)
  
  ## groups variable
  nSubjectGroupsTot <- table(groups)
  wsubjectGroupsTot <- tapply(weights, groups, sum)
  groups <- groups[names(groups) %in% subjectUsed]
  nSubjectGroupsUse <- table(groups)
  wsubjectGroupsUse <- tapply(weightsUse, groups, sum)
  groupsLevs <- levels(groups)
  nGroups <- nlevels(groups)
  groups <- as.integer(groups)
    
  ## order the data
  SUBS <- order(subject, x, y)
  subject <- subject[SUBS]
  x <- x[SUBS]
  y <- y[SUBS]
  
  ## Find curves to plot
    
  ## some needed variables
  x.list <- tapply(X=x,
                   INDEX=list(subject),
                   FUN=function(x){return(x)})  
  trajString <- tapply(X = y, INDEX = list(subject, x), paste,
                        collapse = ",")
  trajString <- apply(X = trajString, MARGIN = 1, FUN = otsplot_data2str)
  
  ## find trajectories
  trajGroup <- as.integer(factor(trajString))
  trajSubject <- 1:nSubjectUse
  groups.name <- 1:max(trajGroup)
  groups.subject <- tapply(trajSubject, trajGroup, function(x){x[1]}) 
  
  ## prepare point data
    
  ## setup pts object
  pts <- data.frame(subject = subject[subject %in% groups.subject],
                    x = x[subject %in% groups.subject],
                    y = y[subject %in% groups.subject])
  pts$traj <- as.integer(factor(pts$subject,
                                levels = groups.subject,
                                labels = groups.name))
  ntraj <- max(pts$traj)
  
  ## if lines should to run from the highest to lowest y category
  ## in case of simultaneous observations
  if (control$plot_gp$lcourse == "downwards") {
    pts <- pts[order(pts$traj, pts$x,
                     factor(pts$y, levels = rev(sort(unique(pts$y))))), ,
               drop = FALSE]
  } 
  
  ## find multiple equal simultaneous observations
  ## and note their frequency
  TMP <- as.data.frame(table(x=pts$x, y=pts$y, traj = pts$traj),
                       responseName = "duplicates")    
  TMP <- TMP[TMP$duplicates != 0,]
  pts <- merge(x = unique(pts), y = TMP, by = c("x", "y", "traj"),
               sort = FALSE)
  
  ## some useful column extractors
  freqCols <- paste("n", 1:nGroups, sep = ".")
  widthCols <- paste("wdt", 1:nGroups, sep = ".")
  colCols <- paste("col", 1:nGroups, sep = ".")
   
  ## weighted point sizes
  pts[, freqCols] <-
    tapply(weightsUse, list(trajGroup, groups), sum)[pts$traj, ]
  pts[is.na(pts)] <- 0
  
  ## jitter coordinates
    
  ## needed variables
  trajWeightsMax <- apply(X = pts[,freqCols, drop = FALSE], MARGIN = 2,
                          FUN = tapply, list(pts$traj), max)
  trajPropUse <- scale(x = matrix(trajWeightsMax, ncol = nGroups),
                       center = FALSE, scale = wsubjectGroupsUse)
  trajPropTot <- scale(x = matrix(trajWeightsMax, ncol = nGroups),
                       center = FALSE, scale = wsubjectGroupsTot)
  trajPropTot <- rbind(trajPropTot, apply(trajPropTot, 2, function(x) 1 - sum(x)))
  trajPropUseMax <- apply(X = trajPropUse, MARGIN = 1, FUN = max)
  trajord <- order(trajPropUseMax, decreasing = TRUE)
  
  otsplot_jitter <- function(data, maxit) {
    
    ## create initial grid
    ngrid <- ngrid0 <- 10
    sl <- ceiling(ngrid * sqrt(trajPropUseMax)) # size
    ngrid <- ceiling(sqrt(sum(sl ^ 2)))
    grid <- matrix(0, ngrid, ngrid) # generate initial grid
    gridsubs <- 1:(ngrid ^ 2) # identifiers
    data$xpos <- data$ypos <- NA      
    potpos <- c()
    blacklist <- c()
    
    ## the algorithm
    for (i in 1:length(trajord)) {
      subscripts <- data$traj == trajord[i]
      count <- 0
      found <- FALSE
      while ((!found)&(count <= maxit)) {
        if ((i > 1) & (count == 0)) {
          if (sl[trajord[i]] < sl[trajord[i-1]]) {
            potpos <- c()
            blacklist <- c()
          }
        }
        count2 <- 0
        while (length(potpos) == 0) { # no potential position available
          if (count2 > 0) { # enlarge grid
            ngrid <- ngrid + 1
            grid <- rbind(grid, rep(0, ngrid - 1))
            grid <- cbind(grid, rep(0, ngrid))
            gridsubs <- seq(1, ngrid^2, 1)
            ## correct blacklist
            blacklist <-  blacklist + floor(blacklist/ngrid)
          }
          ## determine potential grid positions
          potpos <- (grid[gridsubs] == 0) & # occupied positions
          ((gridsubs %% ngrid) < (ngrid - sl[trajord[i]] + 2)) & # to close top
          (gridsubs < (ngrid * (ngrid - sl[trajord[i]] + 1))) # to close left
          if (sl[trajord[i]] > 1) { # filter for groups with ns > 1 field
            potpos <- potpos & (gridsubs %% ngrid != 0)
          }
          potpos <- which(potpos)
          potpos <- potpos[!potpos %in% blacklist]
          count2 <- count2+1
        }     
        xpotpos <- floor(potpos / ngrid)+1
        ypotpos <- potpos %% ngrid
        ypotpos[ypotpos == 0] <- ngrid
        posind <- sample(1:length(potpos), 1) # random assignement
        pos <- potpos[posind]
        xpos <- xpotpos[posind]
        ypos <- ypotpos[posind]
        xsubs <- seq(xpos, xpos + sl[trajord[i]] - 1, 1)
        ysubs <- seq(ypos, ypos + sl[trajord[i]] - 1, 1)
        if (sum(c(grid[ysubs, xsubs]) != 0) == 0) { # assign position
          data$xpos[subscripts] <- xpos
          data$ypos[subscripts] <- ypos
          grid[ysubs,xsubs] <- trajord[i] # register assignement in matrix
          potpos <- potpos[!potpos %in% which(grid == trajord[i])] 
          found <- TRUE
        } else { # nothing found
          count <- count + 1
          blacklist <- c(blacklist, pos)
          potpos <- potpos[-posind]
        }
        if (count == maxit) {
          TMP1 <- data[data$traj == i,]
          TMP2 <- otsplot_data2str(y = TMP$y1, x = TMP$x1)
          warning(paste(" found no grid position for trajectory: ",
                        TMP2 ,", omit!",sep=""))
        }
      }
    }
    data$xjitter <- (data$xpos - 1L) / ngrid
    data$yjitter <- (data$ypos - 1L) / ngrid
    ## determine central plot coordinates
    
    return(list(jitter = data[,c("xjitter", "yjitter")],
                ngrid0 = ngrid0, ngrid = ngrid))
  }
  TMP <- otsplot_jitter(pts, control$maxit)
  pts[,c("xjitter","yjitter")] <- TMP$jitter
  ngrid0 <- TMP$ngrid0
  ngrid <- TMP$ngrid
  
  ## curve colors
  if (is.null(control$plot_gp$col)) { # colors are not specified by the user
    control$plot_gp$col <- rainbowPalette[rep(seq(1, 8), ceiling(ntraj / 8))][1:ntraj]
  } else {
    control$plot_gp$col <- rep(control$plot_gp$col,length.out = ntraj) 
  }
  control$plot_gp$col <- control$plot_gp$col[order(trajord)]
  if (is.character(control$plot_gp$col)) {
    control$plot_gp$col <- col2rgb(control$plot_gp$col)
    control$plot_gp$col <- rgb(control$plot_gp$col[1,],
                               control$plot_gp$col[2,],
                               control$plot_gp$col[3,],
                               control$plot_gp$alpha * 255, maxColorValue = 255)
  }
  
  pts[, colCols] <- matrix(rep(control$plot_gp$col[pts$traj], nGroups), ncol = nGroups)
  control$plot_gp$noobs <- rep(control$plot_gp$noobs, nGroups)
  
  ## color gradients
  if (!is.null(filter)) {
        
    TMP1 <- matrix(NA, ncol = nGroups, nrow = ntraj + 1)
    for (i in 1:ncol(trajPropTot))
      TMP1[, i] <- filter$method(trajPropTot[,i], level = filter$level)
    TMP3 <- TMP1
    for (i in 1:nrow(TMP1)) {
      for (j in 1:ncol(TMP1)) {
        TMP3[i, j] <- otsplot_color(c(0, TMP1[i, j], max(TMP1[, j])), control$plot_gp$hide, c(control$plot_gp$col, control$plot_gp$noobs[1])[i])[2]
      }
    }
    for (i in 1:nGroups) {
      pts[, colCols[i]] <- TMP3[-nrow(TMP3), i][pts$traj]
    }
    control$plot_gp$noobs <- TMP3[nrow(TMP3), ]
  }
  
  ## xlab and ylab
  if (missing(main)) main <- NULL
  if (missing(xlab)) xlab <- deparse(match.call()$x)
  if (missing(ylab)) ylab <- NULL
  
  ## xlim and ylims
  if (missing(xlim)) {
    xlim <- c(1 - sqrt(control$grid_gp$scale) / 2 - sqrt(max((wsubjectGroupsTot - wsubjectGroupsUse) / wsubjectGroupsTot)) * sqrt(control$grid_gp$scale) * ngrid0 / ngrid * control$plot_gp$cex, nXLevs + sqrt(control$grid_gp$scale) / 2)
    xlim <- xlim + c(-1, 1) * 0.025 * diff(range(xlim))
  }
  if (missing(ylim)) {
    ylim <- c(1 - sqrt(control$grid_gp$scale) / 2 - sqrt(max((wsubjectGroupsTot - wsubjectGroupsUse) / wsubjectGroupsTot)) * sqrt(control$grid_gp$scale) * ngrid0 / ngrid * control$plot_gp$cex, nYLevs + sqrt(control$grid_gp$scale) / 2)
    ylim <- ylim + c(-1, 1) * 0.025 * diff(range(ylim))
  }
  
  ## point widths
  pts[, widthCols] <- data.frame(sqrt(scale(x = pts[, freqCols, drop = FALSE], center = FALSE, scale = wsubjectGroupsUse)) * sqrt(control$grid_gp$scale) * ngrid0 / ngrid)
  
  ## data frame for line segments
  ntraj <- max(pts$traj)
  lns <- NULL
  for (i in unique(pts$traj)) {
    SUBS <- which(pts$traj == i)
    nSUBS <- length(SUBS)
    if (length(SUBS) > 1) {
      TMP <- data.frame(traj = rep(i, nSUBS-1),
                        x0 = pts$x[SUBS[-nSUBS]],
                        y0 = pts$y[SUBS[-nSUBS]],
                        x1 = pts$x[SUBS[-1]],
                        y1 = pts$y[SUBS[-1]])
      lns <- rbind(lns,TMP)
    }
  }
  
  if (!is.null(lns)) {
  
    ## merge coordinates and colors
    lns <- merge(x = lns, y = pts[, c("traj", "x", "y", "xjitter", "yjitter")], by.x = c("traj", "x0", "y0"), by.y = c("traj", "x", "y"), all.x = TRUE, all.y = FALSE, sort = FALSE)
    lns <- merge(x = lns, y = pts[,c("traj","x","y","xjitter","yjitter", freqCols, colCols, widthCols)], by.x = c("traj","x1","y1"), by.y = c("traj","x","y"), all.x = TRUE, all.y = FALSE, sort = FALSE)
    colnames(lns)[6:9] <- c("x0jitter", "y0jitter", "x1jitter", "y1jitter")
  }
  
  if (nGroups > 1) { # plot layout
    if (is.null(control$layout)) {
      optpr <- otsplot_layout(nGroups, nGroups, nGroups, 1)
      nLayCols <- optpr[1L]
      nLayRows <- optpr[2L]
    } else {
      nLayCols <- control$layout[2L]
      nLayRows <- control$layout[1L]
    }
  } else { nLayCols <- 1L; nLayRows <- 1L }
  
  ## coordinates for background rectangles
  backrect <- expand.grid(xgrid = 1:nXLevs, ygrid = 1:nYLevs)

  ## reshape pts and lns
  ## pts <- reshape(pts, varying = list(freq = freqCols, width = widthCols, col = colCols), v.names = c("freq", "width", "col"), timevar = "groups", times = 1:nGroups, idvar = "traj", ids = pts$traj, direction = "long", new.row.names = 1:(ntraj * nrow(pts)))
  ## lns <- reshape(lns, varying = list(freq = freqCols, width = widthCols, col = colCols), v.names = c("freq", "width", "col"), timevar = "groups", times = 1:nGroups, idvar = "traj", ids = lns$traj, direction = "long", , new.row.names = 1:(ntraj * nrow(pts)))
  pts <- reshape(pts, varying = list(freq = freqCols, width = widthCols, col = colCols), v.names = c("freq", "width", "col"), timevar = "groups", times = 1:nGroups, direction = "long")
  lns <- reshape(lns, varying = list(freq = freqCols, width = widthCols, col = colCols), v.names = c("freq", "width", "col"), timevar = "groups", times = 1:nGroups, direction = "long")
  
  ## plot data object
  return(structure(
           list(pts = pts, lns = lns, backrect = backrect, 
                trajPropUse = trajPropUse, trajPropUseMax = trajPropUseMax,
                trajPropTot = trajPropTot, ntraj = ntraj,
                nGroups = nGroups, groupsLevs = groupsLevs,
                nXLevs = nXLevs, nYLevs = nYLevs,
                xLevs = xLevs, yLevs = yLevs,
                main = main, xlab = xlab, ylab = ylab,
                xlim = xlim, ylim = ylim,
                nLayCols = nLayCols, nLayRows = nLayRows,
                ngrid = ngrid, ngrid0 = ngrid0,
                grid_gp = gpar(col = control$grid_gp$col,
                  lwd = control$grid_gp$lwd,
                  fill = control$grid_gp$fill),
                grid_scale = control$grid_gp$scale, 
                strip_gp = control$strip_gp,
                cex = control$plot_gp$cex,
                lwd = control$plot_gp$lwd, lorder = control$plot_gp$lorder,
                col_noobs = control$plot_gp$noobs,
                sf_cex = control$plot_gp$sf_cex,
                sf_cex_leaves = control$plot_gp$sf_cex_leaves,
                newpage = control$newpage, margins = control$margins,
                pop = control$pop), class = "otsplot"))
}

print.otsplot <- function(x, ...) {

  if (x$newpage) grid.newpage()

  axTicks <- pretty(1:x$nXLevs)
  axTicks <- intersect(axTicks, 1:x$nXLevs)
  
  if (x$nGroups == 1L) {
    
    pushViewport(plotViewport(xscale = x$xlim,
                              yscale = x$ylim,
                              default.units = "native",
                              margins = x$margins,
                              name = "otsplot"))
    grid.rect()
    otsplot_panel(x, 1L)
    grid.xaxis(axTicks, x$xLevs[axTicks])
    grid.yaxis(1:x$nYLevs, x$yLevs)
    grid.text(x$xlab, y = unit(-3, "lines"))
    grid.text(x$ylab, x = unit(-3, "lines"), rot = 90)
    grid.text(x$main, y = unit(1, "npc") + unit(2, "lines"),
              gp = gpar(fontface = "bold"))
    
    if (x$pop) popViewport(1L) else upViewport(1L)
    
  } else {
  
    strUnit <- unit(2, "strheight", "A")
    pushViewport(plotViewport(layout = grid.layout(x$nLayRows, x$nLayCols),
                              name = x$name, margins = x$margins))
    
    for (i in 1:x$nGroups) {

      ## cell position
      posr <- floor((i - 1) / x$nLayCols) + 1
      posc <- i - floor((i - 1) / x$nLayCols) * x$nLayCols
      
      pushViewport(viewport(layout.pos.col = posc,
                            layout.pos.row = posr,
                            name = paste("otsplot_cell", i, sep = "_"),
                            layout = grid.layout(2, 1,
                              heights = unit.c(strUnit, unit(1, "npc") - strUnit))))  
      
      ## header
      pushViewport(viewport(layout.pos.row = 1, name = paste("strip", i , sep = ".")))
      grid.rect(gp = x$strip_gp)
      grid.text(x$groupsLevs[i], gp = x$strip_gp)
      upViewport(1L)

      ## plot
      pushViewport(viewport(layout.pos.row = 2L, name = paste("plot", i, sep = "_")))
      grid.rect()
      pushViewport(viewport(xscale = x$xlim,
                            yscale = x$ylim,
                            default.units = "native"))
      
      otsplot_panel(x, i)

      if (posc == 1L) grid.yaxis(1:x$nYLevs, x$yLevs)
      if (i > x$nGroups - x$nLayCols) grid.xaxis(axTicks, x$xLevs[axTicks])
      
      upViewport(3L)
    }
    grid.text(x$xlab, y = unit(-3, "lines"))
    grid.text(x$ylab, x = unit(-3, "lines"), rot = 90)
    grid.text(x$main, y = unit(1, "npc") + unit(2, "lines"),
              gp = gpar(fontface = "bold"))
    
    if (x$pop) popViewport(1L) else upViewport(1L)
  }
}


otsplot_panel <- function(x, groups) {

  xf <- diff(x$xlim) / as.numeric(convertWidth(unit(1, "npc"), "inches"))
  yf <- diff(x$ylim) / as.numeric(convertHeight(unit(1, "npc"), "inches"))
  TMP <- max(xf, yf)
  xf <- xf/TMP; yf <- yf/TMP
  
  ## extract data
  x$pts <- x$pts[x$pts$groups == groups, ]
  x$lns <- x$lns[x$lns$groups == groups, ]
    
  ## temporary adaption of line widths
  x$lns$lwd <- 96 * as.numeric(convertUnit(unit(x$lns$width, "native"), "inches", ifelse(xf < yf, "x", "y"), "dimension")) * min(c(xf, yf)) * x$cex * x$lwd
  
  ## translation zones
  grid.rect(x = unit(x$backrect$xgrid, "native"),
            y = unit(x$backrect$ygrid, "native"),
            width = unit(xf * sqrt(x$grid_scale), "native"),
            height = unit(yf * sqrt(x$grid_scale), "native"),
            gp = x$grid_gp)
  
  for (i in order(x$trajPropUse[,groups],
                  decreasing = (x$lorder == "background"))) {
    
    ## extract subscripts
    SUBSpts <- (x$pts$traj %in% i) & (x$pts$freq > 0L)
    SUBSsun <- SUBSpts & (x$pts$duplicates > 1L)
    if (!is.null(x$lns)) {
      SUBSlns <- (x$lns$traj %in% i) & (x$lns$freq > 0L)
      } else {
        SUBSlns <- FALSE
      }    
    
      ## points
    if (sum(SUBSpts) > 0L) {
      grid.rect(x = unit((x$pts$x + xf * (- 0.5 + x$pts$xjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[SUBSpts], "native"), y = unit((x$pts$y + yf * (- 0.5 + x$pts$yjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[SUBSpts], "native"), width = unit((x$pts$width * xf * x$cex)[SUBSpts], "native"), height = unit((x$pts$width * yf * x$cex)[SUBSpts], "native"), gp = gpar(col = x$pts$col[SUBSpts], fill = x$pts$col[SUBSpts], lwd = 0, lex = 0))
      }

    ## lines
    if (sum(SUBSlns) > 0L) {
      grid.segments(x0 = unit((x$lns$x0 + xf * (- 0.5 + x$lns$x0jitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[SUBSlns], "native"), y0 = unit((x$lns$y0 + yf * (- 0.5 + x$lns$y0jitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[SUBSlns], "native"), x1 = unit((x$lns$x1 + xf * (- 0.5 + x$lns$x0jitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[SUBSlns], "native"), y1 = unit((x$lns$y1 + yf * (- 0.5 + x$lns$y0jitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[SUBSlns], "native"), gp = gpar(lwd = x$lns$lwd[SUBSlns], col = x$lns$col[SUBSlns], lineend = "round"))
    }
    
    ## sunflowers
    if (sum(SUBSsun) > 0L) {
      
      grid.points(x = unit((x$pts$x + xf * (- 0.5 + x$pts$xjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[SUBSsun], "native"), y = unit((x$pts$y + yf * (- 0.5 + x$pts$yjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[SUBSsun], "native"), gp = gpar(cex = x$sf_cex), pch = 16)
      i.multi <- which(SUBSsun) # stolen from sunflowerplot()
      ppin <- par("pin")
      pusr <- par("usr")
      xr <- x$pts$width[SUBSsun] / 2 * x$sf_cex_leaves * xf
      yr <- x$pts$width[SUBSsun] / 2 * x$sf_cex_leaves * yf
      i.rep <- rep.int(i.multi, x$pts$duplicates[SUBSsun])
      z <- numeric()
      for (k in i.multi) {
        z <- c(z, 1:x$pts$duplicates[k])
      }
      deg <- (2 * pi * z)/x$pts$duplicates[i.rep]
      
      grid.segments(x0 = unit((x$pts$x + xf * (- 0.5 + x$pts$xjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[i.rep], "native"), y0 = unit((x$pts$y + yf * (- 0.5 + x$pts$yjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[i.rep], "native"), x1 = unit((x$pts$x + xf * (- 0.5 + x$pts$xjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[i.rep] + xr * sin(deg), "native"), y1 = unit((x$pts$y + yf * (- 0.5 + x$pts$yjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[i.rep] + yr * cos(deg), "native"))
    }
  }
}
  

## create trajectory strings
otsplot_data2str <- function(y, x = 1L:length(y),
                         pa.left = "(", pa.right = ")",
                         pot= "^" ,con = "-") {
  subscripts <- !is.na(y);
  ret <- paste(pa.left, y, pa.right, pot, x, sep = "");
  ret <- ret[subscripts];
  ret <- paste(ret, collapse = con)
  return(ret)
}


## layout optimization function 
otsplot_layoutObjFun <- function(nXLevs, nYLevs, n, c = 0.5) {
  control <- nXLevs * nYLevs - n # number of empty windows
  ## optimization value: function of number of empty windows
  ## and ncol/nrow ratio 
  opt <- control / n + c * abs(1 - min(nXLevs, nYLevs) / max(nXLevs, nYLevs))
  return(c(control,opt))
}

## layout optimizer
otsplot_layout <- function(nXLevs, nYLevs, npanels, c = 1) {
  minvalue1 <- otsplot_layoutObjFun(nXLevs - 1, nYLevs, npanels, c)
  minvalue2 <- otsplot_layoutObjFun(nXLevs, nYLevs - 1, npanels, c)
  while (minvalue1[1] >= 0 | minvalue2[1] >= 0) {
    if ((minvalue1[1] >= 0) & (minvalue2[1] >= 0)) {
      if (minvalue1[2] < minvalue2[2]) {
        nXLevs <- nXLevs - 1 } else { nYLevs <- nYLevs - 1 }
    } else {
      if (minvalue1[2] >= 0) {
        nXLevs <- nXLevs - 1
      } else {
        if (minvalue2[2] >= 0) {
          nYLevs <- nYLevs - 1
        }
      }}
    minvalue1 <- otsplot_layoutObjFun(nXLevs - 1, nYLevs, npanels, c)
    minvalue2 <- otsplot_layoutObjFun(nXLevs, nYLevs - 1, npanels, c)
  }
  return(c(nXLevs, nYLevs))
}


## function to colour lines
otsplot_color <- function(value, col1, col2) {
  mp <- colorRamp(c(col1, col2))
  col <- rgb(mp(value), maxColorValue = 255)
}

## linear colour filter function
linear <- function(x, ...) {
  return((x - min(x)) / diff(range(x)))
}

## minimal frequency for lines to filter
minfreq <- function(x, level = 0.05) {
  return(1*(x >= level))
}

## filter a given proportion of most frequent sequences
cumfreq <- function(x, level = 0.75) {
  TMP <- which(cumsum(sort(x, decreasing = TRUE)) >= level)[1]
  ret <- vector("logical", length(x))
  ret[order(x, decreasing = TRUE)[1:TMP]] <- TRUE
  return(1 * ret)
}
