#' meanbranches
#'
#' Calculate the mean branch lengths from a BayesTraits RJ analysis posterior, when topology is fixed.
#' @param trees An object of class multiPhylo - each tree must have the same topology.
#' @param ladderize Ladderize the input trees or not
#' @param returntree Return the mean branch-lengthed tree or not? Defaults to FALSE, and only works if topologies are identical between trees.
#' @export

meanBranches <- function(trees, ladderize = FALSE, returntree = FALSE) {
  # gather all the edge.lengths into a list.
  bls <- list()
  if (ladderize == TRUE) {

    for (i in 1:length(trees)) {
      bls[[i]] <- ladderize(trees[[i]])$edge.length
    }
      
    res <- colMeans(do.call(rbind, bls))
    
    if (returntree) {
      tree <- ladderize(trees[[1]])
      tree$edge.length <- res
      res <- list(bls = bls, tree = tree)
    }
    
  } else {

    for (i in 1:length(trees)) {
      bls[[i]] <- trees[[i]]$edge.length
    }
    
    res <- colMeans(do.call(rbind, bls))
    
    if (returntree) {
      tree <- trees[[1]]
      tree$edge.length <- res
      res <- list(bls = bls, tree = tree)
    }
        
  }
    
  return(res)
}

