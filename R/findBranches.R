#' findBranches
#' Finds the branch numbers that descend from a given node of a tree.
#' @param tree An object of class phylo
#' @param node The node the descendent branches of which are of interest
#' @param tail Include the branch leading to the node of interest, or not?
#' @name findBranches
#' @export

findBranches <- function(tree, node, tail = TRUE) {
  allbranches <- NULL
  for (i in node) {
    descs <- getDescs(tree, i)
    tips <- descs[which(descs <= length(tree$tip.label))]
    internal <- descs[which(descs > length(tree$tip.label))]
    if (tail) {
      allbranches <- c(allbranches, which((tree$edge[ , 2] == i)), which.edge(tree, tree$tip.label[tips]))
    } else {
      allbranches <- c(allbranches, which((tree$edge[ , 1] == i)), which.edge(tree, tree$tip.label[tips]))
      allbranches <- unique(allbranches)
    }
  }
  return(allbranches)
}
