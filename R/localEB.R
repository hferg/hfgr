#' localEB
#'
#' A function to transform a tree according to the exponential change/Early Burst model of Harmon et al 2010
#' @param tree A tree object of class phylo
#' @param node The node that the transformation should be applied at (includes all branches and nodes downstream of here).
#' @param a The EB parameter for the transformation.
#' @param rescale Rescale the tree after transformation? Defailts to TRUE.

localEB <- function(tree, node, a, rescale = TRUE) {
  descs <- getDescs(tree, node)
  descs <- c(descs, node)
  trans.edges <- which(tree$edge[ ,1] %in% descs & tree$edge[ ,1] > length(tree$tip.label))
  bls <- tree$edge.length[trans.edges]
  times <- branching.times(tree)
  times <- times[names(times) %in% descs]
  maxt <- max(times)
  originalsum <- sum(tree$edge.length[trans.edges])
  res <- tree
  for (i in trans.edges) {
    branch <- res$edge.length[i]
    age <- times[which(names(times) == res$edge[i, 1])]
    t1 <- max(times) - age
    t2 <- t1 + branch
    res$edge.length[i] <- (exp(a * t2) - exp(a * t1))/(a)
  }
  
  if (rescale) {
    newsum <- sum(res$edge.length[trans.edges])
    ratio <- originalsum / newsum
    res$edge.length[trans.edges] <- res$edge.length[trans.edges] * ratio
  }
  
  return(res)
}

