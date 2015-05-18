#' localDelta
#'
#' A function that transforms a tree according to lambda, or a portion of a tree according to lambda.
#' @param tree A tree of class phylo.
#' @param node A node number describing the clade(s) to be transformed.
#' @param delta The value or values of lambda by which to transform the specified clade(s).

localDelta <- function(tree, node, delta) {
  descs <- getDescs(tree, node)
  descs <- c(descs, node)
  trans.edges <- which(tree$edge[ ,1] %in% descs & tree$edge[ ,1] > length(tree$tip.label))
  bls <- tree$edge.length[trans.edges]
  times <- branching.times(tree)
  times <- times[names(times) %in% descs]
  original.length <- max(times)
  times <- max(times) - times
  res <- tree
  
  for (i in 1:length(trans.edges)) {
    bl <- tree$edge.length[trans.edges[i]]
    age <- age <- times[names(times) %in% tree$edge[trans.edges][i]]
    res$edge.length[trans.edges[i]] <- (age + bl) ^ delta - age ^ delta
  }
  
  return(res)
}

