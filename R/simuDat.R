#' simuDat
#'
#' Simulate data over a phylogeny, or sub section of a phylogeny, in order to reflect variable rates.
#' Based on adding a constant to the tip data branch by branch and weighting by branch length and does
#' not need to stretch the phylogeny.

#' @param tree
#' @param node
#' @param a The constant that is weighted by branch length and added to the data at the tips. Once the relationship between scalar and this is known, this can be scalar
#' @param dat If you want to use some data you already have and change it this way, put it here.
#' @export

simuDat <- function(tree, node, a, sig, dat = NULL) {
  
  if (isDefined(dat)) {
    dat <- dat
  } else {
    dat <- rTraitCont(tree, model = "BM", sigma = sig)
  }
  
  transdat <- dat  

  descs <- c(getDescs(tree, node), node)
  internal <- descs[descs > length(tree$tip.label)]
  tips <- descs[descs < length(tree$tip.label)]
  edges <- tree$edge[tree$edge[ ,2] %in% c(internal, tips), ]
  bls <- tree$edge.length[tree$edge[ ,2] %in% c(internal, tips)]
  
  for (i in 1:nrow(edges)) {
    tps <- getDescs(tree, edges[i, 2])
    tps <- tps[tps < length(tree$tip.label)]
    num <- sample(2, 1)
    
    if (num == 1) {
      transdat[tps] <- transdat[tps] + (a * (sqrt(sig)))
    } else {
      transdat[tps] <- transdat[tps] - (a * (sqrt(sig)))
    }
  }
  
  dat <- data.frame(taxon = names(dat), x = dat)
  transdat <- data.frame(taxon = names(transdat), x = transdat)
  res <- list(original = dat, transformed = transdat)
  return(res)
}

