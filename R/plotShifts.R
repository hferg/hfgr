#' plotShifts
#' 
#' Plots the locations of the origins of scalars from the postprocessor output of bayestraits.
#' CURRENTLY WORKS ONLY FOR DELTAS.
#' @param PP The psotprocessor (localscalrPP) output.
#' @param scalar The scalar to find and plot from the post processor - delta/lambda/kappa/node/branch
#' @param scaled Plot the original tree (scaled = "time", the default), or the mean/sclaed tree (scaled = "mean") or plot the tree scaled only by scalars present above the threshold (scaled = "threshold")?
#' @param measure When plotting "siginficant" tree, what measure of the parameter? Median (default), mode or mean.
#' @param colour The colour to use for the node circles
#' @param cex The scaling factor for the size of the node circles
#' @param tips Show tip labels?
#' @param threshold Threshold of probability in posterior to display deltas for, defaults to zero (i.e. shows all deltas shaded proportionally to the posterior probability)
#' @param edge.cols A vector of edge colours for the phylogeny
#' @param main Title for the plot
#' @param scale Include scale bar?
#' @param bordercol If using directional arrows, what colour is the border of the triangles?
#' @param direction If TRUE nodes are identified with an upward triangle if the parameter is > 1 and a downward if < 1
#' @param border.width Width of the border of the triangles if using direction.
#' @param exludeones If plotting according to a threshold of significance, should 1s (i.e. no scalar) be excluded from the posterior when calculating average scalar?
#' @param relativetrans If TRUE (defaults to FALSE) the scale of transparency will go from the threshold (totally transparent) to the maximum presence (full opacity).
#' @param transparency Plot the node labels according to their presnce in the posterior?
#' @param gradientcols A vector of two colours - the min and max colours used when colouring the tree according to percentage time rate scaled (when threshold = 0)
#' @name plotShits
#' @export

plotShifts <- function(PP, scalar, threshold = 0, colour = "red", direction = FALSE, 
  scaled = "time",  cex = 1, tips = FALSE, edge.cols = "black", edge.width = 1, main = "", 
  scale = TRUE, bordercol = "black", border.width = 1, measure = "median", excludeones = FALSE,
  relativetrans = FALSE, transparency = TRUE, gradientcols = c("dodgerblue", "firebrick1")) {

  if (scalar == "delta") {
    cl <- "nOrgnDelta"
    par <- paste0(measure, "Delta")
    .mode <- "trans"
  } else if (scalar == "kappa") {
    cl <- "nOrgnKappa"
    par <- paste0(measure, "Kappa")
    .mode <- "trans"
  } else if (scalar == "lambda") {
    cl <- "nOrgnLambda"
    par <- paste0(measure, "Lambda")
    .mode <- "trans"
  } else if (scalar == "rate") {
    cl <- "nOrgnScalar"
    par <- paste0(measure, "Rate")
    mode <- "rate"
  }

  if (mode == "trans") {
    if (threshold != 0) {
      nodes <- PP$data$descNode[which((PP$data[ , cl] / PP$niter) >= threshold)]
      
      if (transparency) {
        alphas <- PP$data[which((PP$data[ , cl] / PP$niter) >= threshold) , cl] / PP$niter
      } else {
        alphas <- rep(1, length(nodes))
      }

      if (relativetrans) {
        for (i in 1:length(alphas)) {
          alphas[i] <- (alphas[i] - min(alphas)) / (max(alphas) - min(alphas))
        }
      }

    } else {
      nodes <- PP$data$descNode[which(PP$data[ , cl] != threshold)]
      
      if (transparency) {
        alphas <- PP$data[which((PP$data[ , cl] / PP$niter) >= threshold) , cl] / PP$niter
      } else {
        alphas <- rep(1, length(nodes))
      }

      if (relativetrans) {
        for (i in 1:length(alphas)) {
          alphas[i] <- (alphas[i] - min(alphas)) / (max(alphas) - min(alphas))
        }
      }

    }

    col <- vector(mode = "character", length = length(nodes))  

    for (j in 1:length(alphas)) {
      col[j] <- makeTransparent(colour, alpha = alphas[j])
    }

    if (direction) {
      shp <- vector(mode = "numeric", length = length(nodes))
      for (i in 1:length(nodes)) {
        if (PP$data[PP$data$descNode == nodes[i], par] > 1) {
          shp[i] <- 24
        } else if (PP$data[PP$data$descNode == nodes[i], par] < 1) {
          shp[i] <- 25
        } else if (PP$data[PP$data$descNode == nodes[i], par] == 1) {
          shp[i] <- 16
        }
      }
    } else {
      shp <- 16
    }

  } else if (mode == "rate") {

    percscaled <- apply(PP$scalars[[1]][2:nrow(PP$scalars[[1]]), ], 1, function(x) sum(x != 1)) / PP$niter

    if (threshold == 0) {
      # Work out the colour ramp.
      # First turn the number of times scaled into percentages.
      edge.cols <- color.scale(percscaled, extremes = gradientcols, na.color = NA)

    } else if (threshold > 0) {
      nodes <- as.numeric(names(percscaled[percscaled >= threshold]))
      edge.cols <- rep("black", nrow(PP$meantree$edge))
      edge.cols[PP$meantree$edge[ , 2] %in% nodes] <- colour
    }

  }

  if (scaled == "time") {
    tree <- PP$meantree
    tree$edge.length <- PP$data$orgBL[2:nrow(PP$data)]
  } else if (scaled == "mean") {
      tree <- PP$meantree
  } else if (scaled == "median") {
    tree <- PP$meantree
    tree$edge.length <- PP$data$medianBL[2:nrow(PP$data)]
  } else if (scaled == "mode") {
    tree <- PP$meantree
    tree$edge.length <- PP$data$modeBL[2:nrow(PP$data)]
  } else if (scaled == "threshold") {
    tree <- significantTransformation(PP = PP, scalar = scalar, threshold = threshold, 
      measure = measure, excludeones = excludeones)
  }

  plotPhylo(tree, tips = tips, edge.cols = edge.cols, edge.width = edge.width, 
    main = main, scale = scale)
  
  if (scalar != "rate") {
    if (direction) {
      nodelabels(node = nodes, bg = col, col = bordercol, pch = shp, 
        cex = cex, lwd = border.width)      
    } else {
      nodelabels(node = nodes, bg = col, pch = shp, cex = cex)  
    }
  }
}
