#' localRate
#'
#' A function that transforms branch lengths according to a scalar.
#' @param tree A tree of class phylo.
#' @param node A node number describing the clade to be transformed.
#' @param scalar The multiplier that rate is increased by in the clade specified by node.

localRate <- function(tree, node, scalar) {
    descs <- getDescs(tree, node)
    descs <- c(node, descs)
    trans.edges <- which(tree$edge[ ,2] %in% descs)
    tree$edge.length[trans.edges] <- scalar * tree$edge.length[trans.edges]
    return(tree)
  }
  
