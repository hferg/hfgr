#' findBranches
#' Finds the branch numbers that descend from a given node of a tree.
#' @param tree An object of class phylo
#' @param node The node the descendent branches of which are of interest
#' @param tail Include the branch leading to the node of interest, or not?
#' @name findBranches
#' @export

findBranches <- function(tree, node, tail = TRUE) {
  descs <- getDescs(tree, node)
  tips <- descs[which(descs <= length(tree$tip.label))]
  internal <- descs[which(descs > length(tree$tip.label))]
  if (tail) {
    allbranches <- c(which((tree$edge[ , 2] == node)), which.edge(tree, tree$tip.label[tips]))
  } else {
    allbranches <- c(which((tree$edge[ , 1] == node)), which.edge(tree, tree$tip.label[tips]))
  }
  
  return(allbranches)
}