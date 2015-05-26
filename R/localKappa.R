#' localKappa
#'
#' A function that transforms branch lengths according to Kappa.
#' @param tree A tree of class phylo.
#' @param node A node number describing the clade to be transformed.
#' @param kappa The multiplier that rate is increased by in the clade specified by node.

# TODO(hfg): Add rescale option.

localKappa <- function(tree, node, kappa, rescale = TRUE) {
  descs <- getDescs(tree, node)
  descs <- c(node, descs)
  trans.edges <- which(tree$edge[ ,2] %in% descs)
  originalsum <- sum(tree$edge.length[trans.edges])
  res <- tree
  res$edge.length[trans.edges] <- tree$edge.length[trans.edges] ^ kappa
  
  if (rescale) {
    newsum <- sum(res$edge.length[trans.edges])
    ratio <- originalsum / newsum
    res$edge.length[trans.edges] <- res$edge.length[trans.edges] * ratio
  }
  
  return(res)
}

