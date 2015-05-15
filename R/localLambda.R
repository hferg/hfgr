#' localLambda
#'
#' A function that transforms a tree according to lambda, or a portion of a tree according to lambda.
#' @param tree A tree of class phylo.
#' @param node A node number or vector of nodes describing the clade(s) to be transformed.
#' @param tips A vector, or list of vectors, of tip labels definining clade(s) to be transformed
#' @param lambda The value or values of lambda by which to transform the specified clade(s).

localLambda <- function(tree, node = NULL, tips = NULL, lambda) {
  lambdaTrans <- function(tree, node, lambda) {
    descs <- getDescs(tree, node)
    trans.edges <- which(tree$edge[ ,2] %in% descs & tree$edge[ ,2] > length(tree$tip.label))
    ht1 <- nodeHeights(tree)
    tree$edge.length[trans.edges] <- lambda * tree$edge.length[trans.edges]
    ht2 <- nodeHeights(tree)
    tree$edge.length[-trans.edges] <- tree$edge.length[-trans.edges] + ht1[-trans.edges, 2] - ht2[-trans.edges, 2]
    return(tree)
  }
  
  if (is.null(tips)) {
   
    if (length(node) != length(lambda)) {
      stop("node and lambda must be the same length")
    }
    
    for (i in 1:length(node)) {
      tree <- lambdaTrans(tree, node[i], lambda[i])
    }
  } else {
  
    if (length(tips) != length(lambda)) {
      stop("tips and lambda must be the same length")
    }
    
    for (i in 1:length(tips)) {
      node <- getMRCA(tree, tips[[i]])
      tree <- lambdaTrans(tree, node, lambda[i])
    }

  }
  
  return(tree)
}

