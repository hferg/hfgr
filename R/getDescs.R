#' getDesc
#'
#' A function to get all descendant nodes from a given node, or vector of tip labels.
#' @param tree A tree of class phylo
#' @param node Either a single node, or a vector of tip labels

getDescs <- function(tree, node, nds = NULL) {

  if (length(node) > 1) {
    node <- getMRCA(tree, node)
  }
  
  if (is.null(nds)) {
    nds <- vector()
  }
    
  dtrs <- tree$edge[which(tree$edge[ , 1] == node), 2]
  nds <- c(nds, dtrs)
  now <- which(dtrs >= length(tree$tip))
  
  if (length(now) > 0) {
  
    for (i in 1:length(now)) {
      nds <- getDescs(tree, dtrs[now[i]], nds)
    }
    
  }
  return(nds)
}

