#' simuDatLambda
#'
#' Simulate data over a phylogeny, or part of a phylogeny, according to the model of Munkemuller et al 2012 by simulating data with a Brownian Motion process, and then adding noise to reflect loss of signal.
#' @param tree A tree of class 'phylo'
#' @param A node that describes the section of the tree to be transformed (if whole tree, the root node).
#' @param sig The standard deviation for the Brownian Motion process.
#' @param w The weight parameter. 1 = strong signal, 0 = weak signal.
#' @param dat A specified dataset, if you want to simulate the loss of signal on an existing dataset, ot something.
#' @export

simuDatLambda <- function(tree, node, sig, w, dat = NULL) {
  
  if (isDefined(dat)) {
    dat <- dat
  } else {
    dat <- rTraitCont(tree, model = "BM", sigma = sig)
  }
  
  transdat <- dat  

  
  descs <- c(getDescs(tree, node), node)
  tips <- descs[descs < length(tree$tip.label)]
  tips <- tree$tip.label[tips]
  
  dd <- transdat[names(transdat) %in% tips]
  dd <- w * dd + (1 - w) * sample(dd)

  transdat[names(transdat) %in% tips] <- dd
  
  res <- list(original = dat, transformed = transdat)
  return(res)
}

