#' plotShifts
#' 
#' Plots the locations of the origins of scalars from the postprocessor output of bayestraits.
#' CURRENTLY WORKS ONLY FOR DELTAS.
#' @param PP The psotprocessor (localscalrPP) output.
#' @param scalar The scalar to find and plot from the post processor - delta/lambda/kappa/node/branch
#' @param scaled Logical - plot the original tree, or the mean/sclaed tree?
#' @param colour The colour to use for the node circles
#' @param cex The scaling factor for the size of the node circles
#' @param tips Show tip labels?
#' @param threshold Threshold of probability in posterior to display deltas for, defaults to zero (i.e. shows all deltas shaded proportionally to the posterior probability)
#' @name plotShits
#' @export

plotShifts <- function(PP, scalar, scaled = FALSE, colour = "black", cex = 1, tips = FALSE, 
  threshold = 0, direction = TRUE, edge.cols = "black") {

  tree <- PP$meantree

  if (scaled == FALSE) {
    tree$edge.length <- PP$data$orgBL[2:nrow(PP$data)]
  }

  if (scalar == "delta") {
    cl <- "nOrgnDelta"
  } else if (scalar == "kappa") {
    cl <- "nOrgnKappa"
  } else if (scalar == "lambda") {
    cl <- "nOrgnLambda"
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

  plotPhylo(tree, tips = tips, edge.cols = edge.cols)
  nodelabels(node = nodes, bg = col, pch = shp, cex = cex)
}
