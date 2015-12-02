#' localDelta
#'
#' A function that transforms a tree according to lambda, or a portion of a tree according to lambda.
#' @param tree A tree of class phylo.
#' @param node A node number describing the clade(s) to be transformed.
#' @param delta The value or values of lambda by which to transform the specified clade(s).
#' @name localDelta

# TODO(hfg): Add rescale option. Find the original root to tip length of the clade, divide it by the new root to tip length, and then multiply the new branch lengths by that.

localDelta <- function(tree, node, delta, rescale = TRUE) {
  descs <- getDescs(tree, node)
  tips <- tree$tip.label[descs[descs <= length(tree$tip.label)]]

  # Make a subtree by using drop.tip, and keep a copy of the edge matrix 
  # for comparison later on...

  subtree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% tips])

  og.edge <- tree$edge[tree$edge[ , 2] %in% descs, ]
  og.edge.length <- tree$edge.length[tree$edge[ , 2] %in% descs]
  
  # And now apply the delta transformation to that new subtree, which means I don't need to 
  # work out transformed edges etc.

  n <- Ntip(subtree)
  hts <- phyloHeights(subtree)
  T <- hts$start[n + 1]
  hts$t <- T - hts$end
  hts$e <- hts$start - hts$end
  hts$a <- hts$t - hts$e

  bls <- (hts$a + hts$e) ^ delta - hts$a ^ delta

  subtree$edge.length <- bls[subtree$edge[ , 2]]

  if (rescale) {
    scale <- T ^ delta
    subtree$edge.length <- (subtree$edge.length / scale) * T
  }

  tree$edge.length[tree$edge[ , 2] %in% descs] <- subtree$edge.length
  res <- tree
  return(res)
}
