#' plotShifts
#' 
#' Plots the locations of the origins of scalars from the postprocessor output of bayestraits.
#' CURRENTLY WORKS ONLY FOR DELTAS.
#' @param PP The psotprocessor (localscalrPP) output.
#' @param scalar The scalar to find and plot from the post processor - delta/lambda/kappa/node/branch
#' @param scaled Plot the original tree (scaled = "time", the default), or the mean/sclaed tree (scaled = "mean") or plot the tree scaled only by scalars present above the threshold (scaled = "threshold")?
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
#' @name plotShits
#' @export

plotShifts <- function(PP, scalar, threshold = 0, colour = "black", direction = TRUE, 
  scaled = "time",  cex = 1, tips = FALSE, edge.cols = "black", edge.width = 1, main = "", 
  scale = TRUE, bordercol = "black", border.width = 1) {

  if (scalar == "delta") {
    cl <- "nOrgnDelta"
    par <- "medianDelta"
    trpar <- "delta"
  } else if (scalar == "kappa") {
    cl <- "nOrgnKappa"
    par <- "medianKappa"
    trpar <- "kappa"
  } else if (scalar == "lambda") {
    cl <- "nOrgnLambda"
    par <- "medianLambda"
    trpar <- "lambda"
  }

  if (threshold != 0) {
    nodes <- PP$data$descNode[which((PP$data[ , cl] / PP$niter) >= threshold)]
    alphas <- PP$data[which((PP$data[ , cl] / PP$niter) >= threshold) , cl] / PP$niter
    values <- PP$data[PP$data$descNode %in% nodes, ]
  } else {
    nodes <- PP$data$descNode[which(PP$data[ , cl] != threshold)]
    alphas <- (PP$data[which(PP$data[ , cl] != 0), cl] / PP$niter)
  }

  col <- vector(mode = "character", length = length(nodes))  

  for (j in 1:length(alphas)) {
    col[j] <- makeTransparent(colour, alpha = alphas[j])
  }

  if (direction) {
    shp <- vector(mode = "numeric", length = length(nodes))
    for (i in 1:length(nodes)) {
      if (PP$data[PP$data$descNode == nodes[i], "medianDelta"] > 1) {
        shp[i] <- 24
      } else if (PP$data[PP$data$descNode == nodes[i], "medianDelta"] < 1) {
        shp[i] <- 25
      }
    }
  } else {
    shp <- 16
  }

  if (scaled == "time") {
    tree <- PP$meantree
    tree$edge.length <- PP$data$orgBL[2:nrow(PP$data)]
  } else if (scaled == "mean") {
      tree <- PP$meantree
  } else if (scaled == "threshold") {
    tree <- PP$meantree
    tree$edge.length <- PP$data$orgBL[2:nrow(PP$data)]    
    dlts <- matrix(nrow = length(nodes), ncol = 3)
    dlts[ , 1] <- nodes
    for (i in 1:nrow(dlts)) {
      dlts[i, 2] <- length(getTipNames(tree, dlts[i, 1]))
      dlts[i, 3] <- PP$data[PP$data$descNode == nodes[i], par]
    }
    dlts <- dlts[order(dlts[ , 2], decreasing = TRUE), ]

    for (i in 1:nrow(dlts)) {
      tree <- treeTrans(tree, node = dlts[i , 1], param = trpar, value = dlts[i, 3])
    }

  }

  plotPhylo(tree, tips = tips, edge.cols = edge.cols, edge.width = edge.width, 
    main = main, scale = scale)
  
  if (direction) {
    nodelabels(node = nodes, bg = col,, col = bordercol, pch = shp, 
      cex = cex, lwd = border.width)      
  } else {
    nodelabels(node = nodes, bg = col, pch = shp, cex = cex)  
  }
}
