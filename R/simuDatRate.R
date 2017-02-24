#' simuDatRate
#'
#' Simulate data over a phylogeny, or sub section of a phylogeny, in order to reflect variable rates.
#' Based on adding a constant to the tip data branch by branch and weighting by branch length and does
#' not need to stretch the phylogeny.
#' @name simuDatRate
#' @param tree The tree that data is simulated on.
#' @param transtree If you want to simulate data using a transformed tree to derive branch lengths, then this is where that tree goes. Defaults to NULL.
#' @param node The node at which to make the changes.
#' @param a The constant that is weighted by branch length and added to the data at the tips. Once the relationship between scalar and this is known, this can be scalar
#' @param sig The standard deviation of the random component of the BM simulation. NOT sigma squared.
#' @param dat If you want to use some data you already have and change it this way, put it here.
#' @param method Either add values recursively by tail branch length (default), or alternative "path" - which adds values proportional to total path length.
#' @param standardise Logical, whether or not to z-standardise the data. Defaults to FALSE.
#' @export

simuDatRate <- function(tree, transtree = NULL, node, a, sig, dat = NULL, method = "recursive", standardise = FALSE) {
  
  if (isDefined(dat)) {
    dat <- dat
  } else {
    dat <- rTraitCont(tree, model = "BM", sigma = sig)
  }
  
  if (standardise == TRUE) {
    dat <- scale(dat, center = TRUE, scale = TRUE)
  }
  
  transdat <- dat

  if (is.null(transtree)) {
    transtree <- tree
  }

  if (node == "root") {
    node <- length(tree$tip.label) + 1
  }

  descs <- c(getDescs(transtree, node), node)
  internal <- descs[descs > length(transtree$tip.label)]
  tips <- descs[descs < length(transtree$tip.label)]
  edges <- transtree$edge[transtree$edge[ ,2] %in% c(internal, tips), ]
  bls <- transtree$edge.length[transtree$edge[ ,2] %in% c(internal, tips)]
  
  if (method == "recursive") {
    
    for (i in 1:nrow(edges)) {
      tps <- getDescs(transtree, edges[i, 2])
      tps <- tps[tps < length(transtree$tip.label)]
      num <- sample(2, 1)
      
      if (num == 1) {
        transdat[tps] <- transdat[tps] + (bls[i] * (a * sig))
      } else {
        transdat[tps] <- transdat[tps] - (bls[i] * (a * sig))
      }
    }
    
  } else if (method == "path") {
  
    print("Not implemented")
  
  }

  res <- list(original = dat, transformed = transdat)
  return(res)
}

