#' meanbranches#'
#' Calculate the mean branch lengths from a BayesTraits RJ analysis posterior, when topology is fixed.
#' @param reftree A tree that provides the reference topology (ideally the tree the analysis was run on)
#' @param trees The posterior sample of stretched trees from which you want the mean branch lengths.
#' @export

meanBranches <- function(reftree, trees, burnin = 0, thinning = 1) {

  reftree <- ladderize(reftree)
  
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
  branches <- vector(mode = "numeric", length(nrow(bls)))
  for (i in 1:nrow(bls)) {
    branches[i] <- mean(bls[i, ])
  }

  for (i in 1:ncol(bls)) {
    bls[ , i] <- bls [ , i] / reftree$edge.length
  }

  for (i in 1:nrow(bls)) {
    quarts <- sort(bls[i, ])
    quart25[i] <- quarts[round(length(quarts) * 0.25)]
    quart75[i] <- quarts[round(length(quarts) * 0.75)]
    rangebl[i] <- max(bls[i, ]) - min(bls[i, ])
  }

  restree <- reftree
  restree$edge.length <- branches
  res <- list(ogtree = reftree, meantree = restree, quart25 = quart25, quart75 = quart75, rangescalar = rangebl)
  return(res)
}
