otsplot <- function(x, ...)
  UseMethod("otsplot")

otsplot.formula <- function(formula, subject, data, subset,
                            na.action = na.pass,
                            drop.unused.levels = FALSE,
                            xlab, ylab, ...) {

  ## checks for group variable
  if (is.character(subject)) subject <- data[, subject]
  if (nrow(data) != length(subject))
    stop("subject vector has not the same length as the data vectors")
  names(subject) <- rownames(data)
  
  ## extract formula
  mc <- match.call(expand.dots = FALSE)
  m <- match(c("data", "subset", "na.action", "drop.unused.levels"),
             names(mc), 0)
  mf <- mc[c(1L, m)]
  mf$formula <- Formula(eval.parent(mc$formula))
  
  ## check formula
  if ((!length(mf$formula)[2] %in% c(1L, 2L)) |
      (length(all.vars(formula(mf$formula, lhs = 0, rhs = 1))) != 1L) |
      (length(all.vars(formula(mf$formula, lhs = 0))) > 2L)) {
    stop("invalid formula")
  }

  ## evaluate model frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval.parent(mf)
  subject <- factor(subject[rownames(mf)])
  
  ## add artifial cond variable
  if (ncol(mf) < 3) {
    cond <- factor(rep(1, nrow(mf)))
  } else {
    cond <- mf[, all.vars(formula)[3]]
  }

  if (missing(xlab)) xlab <- all.vars(formula)[2]
  if (missing(ylab)) ylab <- NULL
    
  ## call default function
  otsplot(x = mf[, all.vars(formula)[2]], y = model.response(mf),
            subject = subject, cond = cond,
            weights = model.weights(mf),
            xlab = xlab, ylab = ylab, ...)
}

otsplot.default <- function(x, y, subject, weights = NULL, cond,
                            plot_gp = gpar(), grid_gp = gpar(),
                            strip_gp = gpar(),
                            filter = NULL, maxit = 200, 
                            main, xlab, ylab, xlim, ylim,
                            layout = NULL, margins = c(5.1, 4.1, 4.1, 3.1),
                            pop = TRUE, newpage = TRUE, ...) {

  ## add default arguments
  
  ## plot parameters
  plot_gp_default <- gpar(cex = 1, lwd = 1 / 4, col = NULL, hide = grey(0.8), noobs = "black", lorder = NULL, lcourse = "upwards", alpha = 1, sf_cex = 0.5, sf_cex_leaves = 1)
  
  ## translation zone parameters
  grid_gp_default <- gpar(scale = 1 / 5, lwd = 1 / 2, fill = grey(0.95), col = grey(0.6))

  ## title font for strips
  strip_gp_default <- gpar(fontsize = 12, fill = grey(0.9))
  
  plot_gp <- appendDefArgs(plot_gp, plot_gp_default)
  grid_gp <- appendDefArgs(grid_gp, grid_gp_default)
  strip_gp <- appendDefArgs(strip_gp, strip_gp_default)
  
  ## check arguments
  
  if (!(plot_gp$lcourse %in% c("upwards", "downwards")))
    stop("invalid lcourse input")
  if (is.null(plot_gp$lorder))
    plot_gp$lorder <- ifelse(is.null(filter), "background", "foreground")    
  if (!(plot_gp$lorder %in% c("background", "foreground")))
    stop("invalid lorder input")
  if (!is.null(filter))
    filter <- otsplot_filter(x = filter)
  
  ## check and prepare raw data
  
  if (sum(c(length(x), length(y)) != length(subject)) != 0) {
    stop("input vectors x, y and subject have different lengths")
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
        warning("problems with distances between x positions. The x positions will not be illustrated adequately.")
        x <- factor(x)
      }
    } else { x <- factor(x) }  # pseudo case
  }
  
  ## y must be an ordered factor
  if (!is.ordered(y)) stop("response must be an ordered factor") 
  
  ## convert subject variable to factor
  if (!is.factor(subject)) subject <- factor(subject)
  
  ## construct/ convert cond variable
  if (!missing(cond)) {
    if (!is.factor(cond)) cond <- factor(cond)
    if (sum(is.na(cond)) > 1) {
      cond <- factor(cond, levels = c(levels(cond), "not available"))
      cond[is.na(cond)] <- "not available"
    }
    if (length(cond) == length(subject)) {
      cond <- unique(data.frame(subject, cond))[, 2]
    }
    if (length(cond) != nlevels(subject)) {
      stop("cannot link cond and subject vector")
    }
    if (nlevels(cond) > length(unique(cond))) {
      warning(paste("erase empty cond levels", paste(setdiff(levels(cond), unique(cond)), collapse = ", ")))
      cond <- factor(cond, levels = levels(cond)[levels(cond) %in% unique(cond)])
    }
  } else cond <- factor(rep(1, nlevels(subject)))
  names(cond) <- levels(subject)
    
  ## construct/ convert weights variable
  if (!is.null(weights)) {
    if (!is.numeric(weights)) stop("weights are not numeric")
    if (min(weights, na.rm = TRUE) < 0) stop("negative weights")
    if (length(weights) == length(subject)) {
      weights <- unique(data.frame(subject, weights))[, 2]
    }
    if (length(weights) != nlevels(subject) | (sum(is.na(weights))) > 0) {
      stop(paste("cannot link weight and subject vector"))
    }
  } else weights <- rep(1, nlevels(subject))

  ## construct/ convert filter variable
  if (!is.null(filter)) {
    if (filter$type == "value") {
      if (!is.numeric(filter$value)) stop("filter$value is not numeric") 
      if (min(filter$value, na.rm = TRUE) < 0) stop("negative filter$value values")
      if (length(filter$value) == length(subject)) {
        filter$value <- unique(data.frame(subject, filter$value))[, 2]
      }
      if (length(filter$value) != nlevels(subject) | (sum(is.na(filter$value))) > 0) {
        stop(paste("cannot link filter$value and subject vector"))
      }
    }
  }

  ## remove NAs and omit.levels in x, y or subject
  SUBS <- is.na(subject) | is.na(x) | is.na(y)
  if (sum(SUBS) > 0) {
    warning(paste(" found NA's in x, y or subject. see entries ", paste(which(SUBS), collapse = ", ")))
    subject <- subject[!SUBS]
    x <- x[!SUBS]
    y <- y[!SUBS]
  }
  
  ## standardize weights
  wsubjectCond.unscaled <- tapply(weights, cond, sum)
  TMP1 <- table(cond)[as.integer(cond)]
  TMP2 <- tapply(weights, cond, sum)[as.integer(cond)]
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
  
  ## cond variable
  nSubjectCondTot <- table(cond)
  wsubjectCondTot <- tapply(weights, cond, sum)
  cond <- cond[names(cond) %in% subjectUsed]
  nSubjectCondUse <- table(cond)
  wsubjectCondUse <- tapply(weightsUse, cond, sum)
  condLevs <- levels(cond)
  nCond <- nlevels(cond)
  cond <- as.integer(cond)
    
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
  cond.name <- 1:max(trajGroup)
  cond.subject <- tapply(trajSubject, trajGroup, function(x){x[1]}) 
  
  ## prepare point data
    
  ## setup pts object
  pts <- data.frame(subject = subject[subject %in% cond.subject],
                    x = x[subject %in% cond.subject],
                    y = y[subject %in% cond.subject])
  pts$traj <- as.integer(factor(pts$subject,
                                levels = cond.subject,
                                labels = cond.name))
  ntraj <- max(pts$traj)
  
  ## if lines should to run from the highest to lowest y category
  ## in case of simultaneous observations
  if (plot_gp$lcourse == "downwards") {
    pts <- pts[order(pts$traj, pts$x, factor(pts$y, levels = rev(sort(unique(pts$y))))), , drop = FALSE]
  } 
  
  ## find multiple equal simultaneous observations
  ## and note their frequency
  TMP <- as.data.frame(table(x=pts$x, y=pts$y, traj = pts$traj), responseName = "duplicates")    
  TMP <- TMP[TMP$duplicates != 0,]
  pts <- merge(x = unique(pts), y = TMP, by = c("x", "y", "traj"),
               sort = FALSE)
  
  ## some useful column extractors
  freqCols <- paste("n", 1:nCond, sep = ".")
  widthCols <- paste("wdt", 1:nCond, sep = ".")
  colCols <- paste("col", 1:nCond, sep = ".")
   
  ## weighted point sizes
  pts[, freqCols] <-
    tapply(weightsUse, list(trajGroup, cond), sum)[pts$traj, ]
  pts[is.na(pts)] <- 0
  
  ## jitter coordinates
    
  ## needed variables
  trajWeightsMax <- apply(X = pts[,freqCols, drop = FALSE], MARGIN = 2,FUN = tapply, list(pts$traj), max)
  trajPropUse <- scale(x = matrix(trajWeightsMax, ncol = nCond), center = FALSE, scale = wsubjectCondUse)
  trajPropTot <- scale(x = matrix(trajWeightsMax, ncol = nCond), center = FALSE, scale = wsubjectCondTot)
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
          if (sl[trajord[i]] > 1) { # filter for cond with ns > 1 field
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
          warning(paste(" found no grid position for trajectory: ", TMP2 ,", omit!",sep=""))
        }
      }
    }
    data$xjitter <- (data$xpos - 1) / ngrid
    data$yjitter <- (data$ypos - 1) / ngrid
    ## determine central plot coordinates
    
    return(list(jitter = data[,c("xjitter", "yjitter")], ngrid0 = ngrid0, ngrid = ngrid))
  }
  TMP <- otsplot_jitter(pts, maxit)
  pts[,c("xjitter","yjitter")] <- TMP$jitter
  ngrid0 <- TMP$ngrid0
  ngrid <- TMP$ngrid
  
  ## curve colors
  if (is.null(plot_gp$col)) { # colors are not specified by the user
    plot_gp$col <- rainbowPalette[rep(seq(1, 8), ceiling(ntraj / 8))][1:ntraj]
  } else {
    plot_gp$col <- rep(plot_gp$col,length.out = ntraj) # one colour for each trajectory type
  }
  plot_gp$col <- plot_gp$col[order(trajord)]
  if (is.character(plot_gp$col)) {
    plot_gp$col <- col2rgb(plot_gp$col)
    plot_gp$col <- rgb(plot_gp$col[1,], plot_gp$col[2,], plot_gp$col[3,], plot_gp$alpha * 255, maxColorValue = 255)
  }
  
  pts[, colCols] <- matrix(rep(plot_gp$col[pts$traj], nCond), ncol = nCond)
  plot_gp$noobs <- rep(plot_gp$noobs, nCond)
  
  ## color gradients
  if (!is.null(filter)) {
        
    if (filter$type == "function") {        
      TMP1 <- matrix(NA, ncol = nCond, nrow = ntraj + 1)
      CALL <- list(filter$value)
      CALL <- append(CALL,filter$args)
      mode(CALL) <- "call"
      for (i in 1:ncol(trajPropTot)) {
        CALL$x <- trajPropTot[,i]
        TMP1[, i] <- eval(CALL)
      }
      TMP3 <- TMP1
      for (i in 1:nrow(TMP1)) {
        for (j in 1:ncol(TMP1)) {
          TMP3[i, j] <- otsplot_color(c(0, TMP1[i, j], max(TMP1[, j])), plot_gp$hide, c(plot_gp$col, plot_gp$noobs[1])[i])[2]
        }
      }
      for (i in 1:nCond) {
        pts[, colCols[i]] <- TMP3[-nrow(TMP3), i][pts$traj]
      }
      plot_gp$noobs <- TMP3[nrow(TMP3), ]
    }
    if (is.vector(filter)) {
      if (is.numeric(filter$value)) {
        
        ## colour by a vector with one numeric value for each individual
        filter$value <- tapply(filter$value[trajSubject], list(trajGroup, cond[trajSubject]), function(x) { median(x) })
        if (is.list(filter$value)) stop("filter vector could not be merged")
        if (nCond == 1) filter$value <- matrix(filter$value, ncol = 1)
        filter$value <- scale(x = filter$value, center = apply(filter$value, 2, min, na.rm = TRUE), scale = apply(filter$value, 2, function(x) { diff(range(x, na.rm = TRUE)) }))
        filter$value[is.na(filter$value)] <- 1
        for (i in 1:ntraj) {
          for (j in 1:nCond) {
            pts[pts$traj == i, colCols[j]] <- otsplot_color(filter$value[i,j], plot_gp$col[i], plot_gp$hide)
          }
        }
      }
    }
  }
  
  ## xlab and ylab
  if (missing(main)) main <- NULL
  if (missing(xlab)) xlab <- match.call()[["x"]]
  if (missing(ylab)) ylab <- NULL
  
  ## xlim and ylims
  if (missing(xlim)) {
    xlim <- c(1 - sqrt(grid_gp$scale) / 2 - sqrt(max((wsubjectCondTot - wsubjectCondUse) / wsubjectCondTot)) * sqrt(grid_gp$scale) * ngrid0 / ngrid * plot_gp$cex, nXLevs + sqrt(grid_gp$scale) / 2)
    xlim <- xlim + c(-1, 1) * 0.025 * diff(range(xlim))
  }
  if (missing(ylim)) {
    ylim <- c(1 - sqrt(grid_gp$scale) / 2 - sqrt(max((wsubjectCondTot - wsubjectCondUse) / wsubjectCondTot)) * sqrt(grid_gp$scale) * ngrid0 / ngrid * plot_gp$cex, nYLevs + sqrt(grid_gp$scale) / 2)
    ylim <- ylim + c(-1, 1) * 0.025 * diff(range(ylim))
  }
  
  ## point widths
  pts[, widthCols] <- data.frame(sqrt(scale(x = pts[, freqCols, drop = FALSE], center = FALSE, scale = wsubjectCondUse)) * sqrt(grid_gp$scale) * ngrid0 / ngrid)
  
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
  
  if (nCond > 1) { # plot layout
    if (is.null(layout)) {
      optpr <- otsplot_layout(nCond, nCond, nCond, 1)
      nLayCols <- optpr[1L]
      nLayRows <- optpr[2L]
    } else {
      nLayCols <- layout[2L]
      nLayRows <- layout[1L]
    }
  } else { nLayCols <- 1L; nLayRows <- 1L }
  
  ## coordinates for background rectangles
  backrect <- expand.grid(xgrid = 1:nXLevs, ygrid = 1:nYLevs)

  ## reshape pts and lns
  pts <- reshape(pts, varying = list(freq = freqCols, width = widthCols, col = colCols), v.names = c("freq", "width", "col"), timevar = "cond", times = 1:nCond, idvar = "traj", ids = pts$traj, direction = "long", new.row.names = 1:(ntraj * nrow(pts)))
  lns <- reshape(lns, varying = list(freq = freqCols, width = widthCols, col = colCols), v.names = c("freq", "width", "col"), timevar = "cond", times = 1:nCond, idvar = "traj", ids = lns$traj, direction = "long", , new.row.names = 1:(ntraj * nrow(pts)))
  
  ## plot data object
  x <- list(pts = pts, lns = lns, backrect = backrect, 
            trajPropUse = trajPropUse, trajPropUseMax = trajPropUseMax,
            trajPropTot = trajPropTot, ntraj = ntraj,
            nCond = nCond, condLevs = condLevs,
            nXLevs = nXLevs, nYLevs = nYLevs,
            xLevs = xLevs, yLevs = yLevs,
            main = main, xlab = xlab, ylab = ylab,
            xlim = xlim, ylim = ylim,
            nLayCols = nLayCols, nLayRows = nLayRows,
            ngrid = ngrid, ngrid0 = ngrid0,
            grid_gp = gpar(col = grid_gp$col, lwd = grid_gp$lwd,
              fill = grid_gp$fill),
            grid_scale = grid_gp$scale, 
            strip_gp = strip_gp,
            cex = plot_gp$cex,
            lwd = plot_gp$lwd, lorder = plot_gp$lorder,
            col_noobs = plot_gp$noobs,
            sf_cex = plot_gp$sf_cex, sf_cex_leaves = plot_gp$sf_cex_leaves,
            newpage = newpage, margins = margins,
            pop = pop)
  class(x) <- "otsplot"
  x
}

print.otsplot <- function(x, ...) {

  if (x$newpage) grid.newpage()

  axTicks <- pretty(1:x$nXLevs)
  axTicks <- intersect(axTicks, 1:x$nXLevs)
  
  if (x$nCond == 1L) {
    
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
    grid.text(x$main, y = unit(1, "npc") + unit(2, "lines"), gp = gpar(fontface = "bold"))
    
    if (x$pop) popViewport(1L) else upViewport(1L)
    
  } else {
  
    strUnit <- unit(2, "strheight", "A")
    pushViewport(plotViewport(layout = grid.layout(x$nLayRows, x$nLayCols), name = x$name, margins = x$margins))
    
    for (i in 1:x$nCond) {

      ## cell position
      posr <- floor((i - 1) / x$nLayCols) + 1
      posc <- i - floor((i - 1) / x$nLayCols) * x$nLayCols
      
      pushViewport(viewport(layout.pos.col = posc, layout.pos.row = posr, name = paste("otsplot_cell", i, sep = "_"), layout = grid.layout(2, 1, heights = unit.c(strUnit, unit(1, "npc") - strUnit))))  
      
      ## header
      pushViewport(viewport(layout.pos.row = 1, name = paste("strip", i , sep = ".")))
      grid.rect(gp = x$strip_gp)
      grid.text(x$condLevs[i], gp = x$strip_gp)
      upViewport(1L)

      ## plot
      pushViewport(viewport(layout.pos.row = 2, name = paste("plot", i, sep = "_")))
      grid.rect()
      pushViewport(viewport(xscale = x$xlim,
                            yscale = x$ylim,
                            default.units = "native"))
      
      otsplot_panel(x, i, ...)

      if (posc == 1L) grid.yaxis(1:x$nYLevs, x$yLevs)
      if (i > x$nCond - x$nLayCols) grid.xaxis(axTicks, x$xLevs[axTicks])
      
      upViewport(3L)
    }
    grid.text(x$xlab, y = unit(-3, "lines"))
    grid.text(x$ylab, x = unit(-3, "lines"), rot = 90)
    grid.text(x$main, y = unit(1, "npc") + unit(2, "lines"), gp = gpar(fontface = "bold"))
    
    if (x$pop) popViewport(1L) else upViewport(1L)
  }
}


otsplot_panel <- function(x, cond) {


  xf <- diff(x$xlim) / as.numeric(convertWidth(unit(1, "npc"), "inches"))
  yf <- diff(x$ylim) / as.numeric(convertHeight(unit(1, "npc"), "inches"))
  TMP <- max(xf, yf)
  xf <- xf/TMP; yf <- yf/TMP
  
  ## extract data
  x$pts <- x$pts[x$pts$cond == cond, ]
  x$lns <- x$lns[x$lns$cond == cond, ]
    
  ## temporary adaption of line widths
  x$lns$lwd <- 96 * as.numeric(convertUnit(unit(x$lns$width, "native"), "inches", ifelse(xf < yf, "x", "y"), "dimension")) * min(c(xf, yf)) * x$cex * x$lwd
  
  ## translation zones
  grid.rect(x = unit(x$backrect$xgrid, "native"),
            y = unit(x$backrect$ygrid, "native"),
            width = unit(xf * sqrt(x$grid_scale), "native"),
            height = unit(yf * sqrt(x$grid_scale), "native"),
            gp = x$grid_gp)
  
  for (i in order(x$trajPropUse[,cond],
                  decreasing = (x$lorder == "background"))) {
    
    ## extract subscripts
    SUBSpts <- (x$pts$traj %in% i) & (x$pts$freq > 0)
    SUBSsun <- SUBSpts & (x$pts$duplicates > 1)
    if (!is.null(x$lns)) {
      SUBSlns <- (x$lns$traj %in% i) & (x$lns$freq > 0)
      } else {
        SUBSlns <- FALSE
      }
    
    
      ## points
    if (sum(SUBSpts) > 0) {
      grid.rect(x = unit((x$pts$x + xf * (- 0.5 + x$pts$xjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[SUBSpts], "native"), y = unit((x$pts$y + yf * (- 0.5 + x$pts$yjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[SUBSpts], "native"), width = unit((x$pts$width * xf * x$cex)[SUBSpts], "native"), height = unit((x$pts$width * yf * x$cex)[SUBSpts], "native"), gp = gpar(col = x$pts$col[SUBSpts], fill = x$pts$col[SUBSpts], lwd = 0, lex = 0))
      }

    ## lines
    if (sum(SUBSlns) > 0) {
      grid.segments(x0 = unit((x$lns$x0 + xf * (- 0.5 + x$lns$x0jitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[SUBSlns], "native"), y0 = unit((x$lns$y0 + yf * (- 0.5 + x$lns$y0jitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[SUBSlns], "native"), x1 = unit((x$lns$x1 + xf * (- 0.5 + x$lns$x0jitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[SUBSlns], "native"), y1 = unit((x$lns$y1 + yf * (- 0.5 + x$lns$y0jitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropUseMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[SUBSlns], "native"), gp = gpar(lwd = x$lns$lwd[SUBSlns], col = x$lns$col[SUBSlns], lineend = "round"))
    }
    
    ## sunflowers
    if (sum(SUBSsun) > 0) {
      
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
otsplot_data2str <- function(y, x = 1:length(y),
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

## convert input for filter argument
otsplot_filter <- function(x) {
  if (is.list(x)) {
    if (sum(names(x) %in% c("type","value")) != 2)
      stop("filter must contain a type and a value object")
    if (!x$type %in% c("value", "function"))
      stop("unknown input for filter$type")
    if ((x$type == "function") & (length(x) >= 2)) {
      if (is.character(x$value)) {
        if (x$value == "linear") {
          x$value <- linear
        } else if (x$value == "minfreq") {
          x$value <- minfreq
        } else if (x$value == "cumfreq") {
          x$value <- cumfreq
        } else {
          stop("invalid x$value argument")
        }
      }
      ret <- list(type = x$type, value = x$value,
                  args = x[!names(x) %in% c("type", "value")])
    } else {
      ret <- x
    }
  }
  if (is.function(x)) {
    ret <- list(type = "function", value = x)
  }
  if (is.vector(x)) {
    if (is.numeric(x)) {
      ret <- list(type = "value", value = x)
    }
    if (is.character(x)) { # could either be a predefined function
      if (x == "linear") {
        ret <- list(type = "function", value = linear)
      } else if (x == "minfreq") {
        ret <- list(type = "function", value = minfreq)
      } else if (x == "cumfreq") {
        ret <- list(type = "function", value = cumfreq)
      }
    }
  }
  return(ret)
}


## function to colour lines
otsplot_color <- function(value, col1, col2) {
  mp <- colorRamp(c(col1, col2))
  col <- rgb(mp(value), maxColorValue = 255)
}

## linear colour filter function
linear <- function(x) {
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
