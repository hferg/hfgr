#' getTipNames
#'
#' A function to get the names of the descendant tips from a given node of a tree.
#' @param tree A tree of class phylo.
#' @param node The node number of interest.
#' export

getTipNames <- function(tree, node) {
  descs <- getDescs(tree, node)
  tree$tip.label[descs <= length(tree$tip.label)]
}
