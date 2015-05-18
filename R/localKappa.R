#' localKappa
#'
#' A function that transforms branch lengths according to Kappa.
#' @param tree A tree of class phylo.
#' @param node A node number of vector of nodes describing the clade to be transformed.
#' @param tips A vector of tip names, or a list of vectors of tip names, that define the areas to be transformed. If used, node is ignored.
#' @param kappa The multiplier that rate is increased by in the clade specified by node.

localKappa <- function(tree, node = NULL, tips = NULL, kappa) {
  rateTrans <- function(tree, node, scalar) {
    descs <- getDescs(tree, node)
    trans.edges <- which(tree$edge[ ,2] %in% descs)
    tree$edge.length[trans.edges] <- tree$edge.length[trans.edges] ^ kappa
    tree$edge.length[which(tree$edge[ ,2] == node)] <- tree$edge.length[which(tree$edge[ ,2] == node)] ^ kappa
    return(tree)
  }
  
  if (is.null(tips)) {
   
    if (length(node) != length(kappa)) {
      stop("node and kappa must be the same length")
    }
    
    for (i in 1:length(node)) {
      
      if (node[i] == length(tree$tip.label) + 1) {
        tree$edge.length <- tree$edge.length ^ kappa[i]
      } else {
        tree <- rateTrans(tree, node[i], kappa[i])
      }
      
    }
  } else {
    
    if (length(tips) != length(kappa)) {
      stop("tips and kappa must be the same length")
    }
    
    for (i in 1:length(tips)) {
      node <- getMRCA(tree, tips[[i]])
      
      if (node == length(tree$tip.label) + 1) {
        tree$edge.length <- tree$edge.length ^ kappa[i]
      } else {
        tree <- rateTrans(tree, node, kappa[i])
      }
      
    }
  }
  
  return(tree)
}

