#' meanbranches
#'
#' Calculate the mean branch lengths from a BayesTraits RJ analysis posterior, when topology is fixed.
#' @param reftree A tree that provides the reference topology (ideally the tree the analysis was run on)
#' @param trees The posterior sample of stretched trees from which you want the mean branch lengths.
#' @export

meanBranches <- function(reftree, trees, burnin = 0, thinning = 1) {

  reftree <- ladderize(reftree)
  
  trees <- trees[seq.int(burnin, length(trees), thinning)]
  
  bls <- vector(mode = "numeric", length = length(reftree$edge.length))
  
  for (i in 1:length(trees)) {
    tree <- ladderize(trees[[i]])
    
    if (sum(reftree$tip.label == tree$tip.label) != length(reftree$tip.label)) {
      stop(paste("Tip labels on tree", i, "do not mactch reference tree"))
    }
    
    if (sum(reftree$edge == tree$edge) != length(reftree$edge)) {
      stop(paste("Tree", i, "has a different topology to reference tree"))
    }
    
    bls <- bls + tree$edge.length
    
  }
  
  bls <- bls / length(trees)
  restree <- reftree
  restree$edge.length <- bls
  res <- list(ogtree = reftree, meantree = restree)
  return(res)
}

