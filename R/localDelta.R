#' localDelta
#'
#' A function that transforms a tree according to lambda, or a portion of a tree according to lambda.
#' @param tree A tree of class phylo.
#' @param node A node number describing the clade(s) to be transformed.
#' @param delta The value or values of lambda by which to transform the specified clade(s).

# TODO(hfg): Add rescale option. Find the original root to tip length of the clade, divide it by the new root to tip length, and then multiply the new branch lengths by that.

localDelta <- function(tree, node, delta, rescale = TRUE) {
  descs <- getDescs(tree, node)
  descs <- c(descs, node)
  trans.edges <- which(tree$edge[ ,1] %in% descs & tree$edge[ ,1] > length(tree$tip.label))
  trans.edges <- c(trans.edges, which(tree$edge[ ,2] == node))
  old.mean <- mean(tree$edge.length[trans.edges])
  descs <- c(descs, tree$edge[which(tree$edge[ ,2] == node), 1])
  bls <- tree$edge.length[trans.edges]
  times <- branching.times(tree)
  times <- times[names(times) %in% descs]
  transnode <- nodeheight(tree, node)
  originalsum <- sum(tree$edge.length[trans.edges])
  times <- max(times) - times
  res <- tree
  
  for (i in 1:length(trans.edges)) {
    bl <- tree$edge.length[trans.edges[i]]
    age <- times[names(times) %in% tree$edge[trans.edges][i]]
    res$edge.length[trans.edges[i]] <- (age + bl) ^ delta - age ^ delta
  }
  
  if (rescale) {
    newsum <- sum(res$edge.length[trans.edges])
    ratio <- originalsum / newsum
    res$edge.length[trans.edges] <- res$edge.length[trans.edges] * ratio
  }
  
  return(res)
}

