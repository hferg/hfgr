#' meanbranches
#'
#' Calculate the mean branch lengths from a BayesTraits RJ analysis posterior, when topology is fixed.
#' @param reftree A tree that provides the reference topology (ideally the tree the analysis was run on)
#' @param trees The logfile for the rj trees posterior.
#' @export

meanBranches <- function(reftree, trees, burnin = 0, thinning = 1) {

  reftree <- ladderize(reftree)
  trees <- read.nexus(trees)
  trees <- trees[seq.int(burnin, length(trees), thinning)]
  
  #bls <- vector(mode = "numeric", length = length(reftree$edge.length))
  bls <- matrix(nrow = length(reftree$edge.length), ncol = length(trees))
  
  for (i in 1:length(trees)) {
    tree <- ladderize(trees[[i]])
    
    if (sum(reftree$tip.label == tree$tip.label) != length(reftree$tip.label)) {
      stop(paste("Tip labels on tree", i, "do not mactch reference tree"))
    }
    
    if (sum(reftree$edge == tree$edge) != length(reftree$edge)) {
      stop(paste("Tree", i, "has a different topology to reference tree"))
    }
    
    bls[ , i] <- tree$edge.length
    
  }
  
  rangebl <- vector(mode = "numeric", length = nrow(bls)) 
  quart25 <- vector(mode = "numeric", length = nrow(bls))
  quart75 <- vector(mode = "numeric", length = nrow(bls))
  branchesmean <- vector(mode = "numeric", length(nrow(bls)))
  branchesmedian <- vector(mode = "numeric", length(nrow(bls)))
  for (i in 1:nrow(bls)) {
    branchesmean[i] <- mean(bls[i, ])
    branchesmedian[i] <- median(bls[i, ])
    quarts <- sort(bls[i, ])
    quart25[i] <- quarts[round(length(quarts) * 0.25)]
    quart75[i] <- quarts[round(length(quarts) * 0.75)]
    rangebl[i] <- max(bls[i, ]) - min(bls[i, ])
  }

  res <- list(ogtree = reftree, 
              meanbranches = branchesmean, 
              medianbranches = branchesmedian, 
              quart25 = quart25, 
              quart75 = quart75, 
              rangescalar = rangebl)
  return(res)
}
