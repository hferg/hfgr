#' significantTransformation
#' 
#' Transforms a tree according to only the scalars present above a certain threshold.
#' @param PP The psotprocessor (localscalrPP) output.
#' @param scalar The scalar to find and plot from the post processor - delta/lambda/kappa/node/branch
#' @param measure Median (default), mean or mode scalar.
#' @param threshold Threshold of probability in posterior to display deltas for, defaults to zero (i.e. shows all deltas shaded proportionally to the posterior probability)
#' @name significantTransformation
#' @export

significantTransformation <- function(PP, scalar, measure = "median", threshold = 0) {

  if (scalar == "delta") {
    cl <- "nOrgnDelta"
    par <- paste0(measure, "Delta")
    trpar <- "delta"
  } else if (scalar == "kappa") {
    cl <- "nOrgnKappa"
    par <- paste0(measure, "Kappa")
    trpar <- "kappa"
  } else if (scalar == "lambda") {
    cl <- "nOrgnLambda"
    par <- paste0(measure, "Lambda")
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

  return(tree)
}
