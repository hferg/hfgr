#' meanbranches
#'
#' Calculate the mean branch lengths from a BayesTraits RJ analysis posterior, when topology is fixed.
#' @param trees An object of class multiPhylo - each tree must have the same topology.
#' @export

meanBranches <- function(trees) {
  # gather all the edge.lengths into a list.
  bls <- list()
  for (i in length(trees)) {
    bls[[i]] <- trees[[i]]$edge.length
  }
  return(colMeans(do.call(rbind, bls)))
}

