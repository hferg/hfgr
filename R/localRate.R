#' localLambda
#'
#' A function that transforms branch lengths according to a scalar.
#' @param tree A tree of class phylo.
#' @param node A node number of vector of tip labels describing the clade to be transformed.
#' @param scalar The multiplier that rate is increased by in the clade specified by node.

localRate <- function(tree, node, scalar) {
  
  if (node == length(tree$tip.label) + 1) {
    tree$edge.length <- tree$edge.length * scalar
  } else {
  descs <- getDescs(tree, node)
  trans.edges <- which(tree$edge[ ,2] %in% descs)
  tree$edge.length[trans.edges] <- scalar * tree$edge.length[trans.edges]
  tree$edge.length[which(tree$edge[ ,2] == node)] <- scalar * tree$edge.length[which(tree$edge[ ,2] == node)]
  }
  
  return(tree)
}



