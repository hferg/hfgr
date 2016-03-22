#' plotShifts
#' 
#' Plots the locations of the origins of scalars from the postprocessor output of bayestraits.
#' CURRENTLY WORKS ONLY FOR DELTAS.
#' @param PP The psotprocessor (localscalrPP) output.
#' @param tree The original phylogeny used for the analysis.
#' @param scalar The scalar to find and plot from the post processor - delta/lambda/kappa/node/branch
#' @param scaled Logical - plot the original tree, or the mean/sclaed tree?
#' @param colour The colour to use for the node circles
#' @param cex The scaling factor for the size of the node circles
#' @param tips Show tip labels?
#' @name plotShits
#' @export

plotShifts <- function(PP, tree, scalar, scaled = FALSE, colour = "black", cex = 1, tips = FALSE) {

  tree <- ladderize(tree)

  if (scaled) {
    tree$edge.length <- PP$data$meanBL
  }

  if (scalar == "delta") {
    cl <- "nOrgnDelta"
  } else if (scalar == "kappa") {
    cl <- "nOrgnKappa"
  } else if (scalar == "lambda") {
    cl <- "nOrgnLambda"
  }

  nodes <- PP$data$descNode[which(PP$data[ , cl] != 0)]
  col <- vector(mode = "character", length = length(nodes))
  alphas <- (PP$data[which(PP$data[ , cl] != 0), "nOrgnDelta"] / PP$niter)

  for (j in 1:length(alphas)) {
    col[j] <- makeTransparent(colour, alpha = alphas[j])
  }

  plotTree(tree, tips = tips)
  nodelabels(node = nodes, col = col, pch = 16, cex = cex)
}
