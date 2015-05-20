#' meanbranches
#'
#' Calculate the mean branch lengths from a BayesTraits RJ analysis posterior, when topology is fixed.
#' @param trees An object of class multiPhylo - each tree must have the same topology.
#' @export

meanBranches <- function(trees, ladderize = FALSE) {
  # gather all the edge.lengths into a list.
  bls <- list()
  if (ladderize == TRUE) {

    for (i in 1:length(trees)) {
      bls[i] <- ladderize(trees[[i]])$edge.length
    }
    
  } else {

    for (i in 1:length(trees)) {
      bls[[i]] <- trees[[i]]$edge.length
    }
      
  }

  return(colMeans(do.call(rbind, bls)))
}

