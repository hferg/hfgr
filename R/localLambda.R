#' localLambda
#'
#' A function that transforms a tree according to lambda, or a portion of a tree according to lambda.
#' @param tree A tree of class phylo.
#' @param node A node number of vector of tip labels describing the clade to be transformed.
#' @param lambda The value of lambda by which to transform the tree.

localLambda <- function(tree, node, lambda) {
  descs <- getDescs(tree, node)
  trans.edges <- which(tree$edge[ ,2] %in% descs & tree$edge[ ,2] > length(tree$tip.label))
  ht1 <- nodeHeights(tree)
  tree$edge.length[trans.edges] <- lambda * tree$edge.length[trans.edges]
  ht2 <- nodeHeights(tree)
  tree$edge.length[-trans.edges] <- tree$edge.length[-trans.edges] + ht1[-trans.edges, 2] - ht2[-trans.edges, 2]
  return(tree)
}

