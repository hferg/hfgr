#' localLambda
#'
#' A function that transforms branch lengths according to a scalar.
#' @param tree A tree of class phylo.
#' @param node A node number of vector of nodes describing the clade to be transformed.
#' @param tips A vector of tip names, or a list of vectors of tip names, that define the areas to be transformed. If used, node is ignored.
#' @param scalar The multiplier that rate is increased by in the clade specified by node.

localRate <- function(tree, node = NULL, tips = NULL, scalar) {
  rateTrans <- function(tree, node, scalar) {
    descs <- getDescs(tree, node)
    trans.edges <- which(tree$edge[ ,2] %in% descs)
    tree$edge.length[trans.edges] <- scalar * tree$edge.length[trans.edges]
    tree$edge.length[which(tree$edge[ ,2] == node)] <- scalar[i] * tree$edge.length[which(tree$edge[ ,2] == node)]
    return(tree)
  }
  
  if (is.null(tips)) {
   
    if (length(node) != length(scalar)) {
      stop("node and scalar must be the same length")
    }
    
    for (i in 1:length(node)) {
      
      if (node[i] == length(tree$tip.label) + 1) {
        tree$edge.length <- tree$edge.length * scalar[i]
      } else {
        tree <- ratetrans(tree, node[i], scalar[i]
      }
      
    }
  } else {
    
    if (length(tips) != length(scalar)) {
      stop("tips and scalar must be the same length")
    }
    
    for (i in 1:length(tips)) {
      node <- getMRCA(tree, tips[[i]])
      
      if (node == length(tree$tip.label) + 1) {
        tree$edge.length <- tree$edge.length * scalar[i]
      } else {
        tree <- ratetrans(tree, node, scalar[i]
      }
      
    }
  }
  
  return(tree)
}

