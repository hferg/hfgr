##############################################################################
#' localRate
#'
#' A function that transforms branch lengths according to a scalar.
#' @param tree A tree of class phylo.
#' @param node A node number describing the clade to be transformed.
#' @param scalar The multiplier that rate is increased by in the clade specified by node.

localRate <- function(tree, node, scalar) {
    descs <- getDescs(tree, node)
    descs <- c(node, descs)
    trans.edges <- which(tree$edge[ ,2] %in% descs)
    tree$edge.length[trans.edges] <- scalar * tree$edge.length[trans.edges]
    return(tree)
  }
  
##############################################################################
#' localLambda
#'
#' A function that transforms a tree according to lambda, or a portion of a tree according to lambda.
#' @param tree A tree of class phylo.
#' @param node A node number describing the clade(s) to be transformed.
#' @param lambda The value or values of lambda by which to transform the specified clade(s).

localLambda <- function(tree, node, lambda) {
    descs <- getDescs(tree, node)
    trans.edges <- which(tree$edge[ ,2] %in% descs & tree$edge[ ,2] > length(tree$tip.label))
    ht1 <- nodeHeights(tree)
    tree$edge.length[trans.edges] <- lambda * tree$edge.length[trans.edges]
    ht2 <- nodeHeights(tree)
    tree$edge.length[-trans.edges] <- tree$edge.length[-trans.edges] + ht1[-trans.edges, 2] - ht2[-trans.edges, 2]
    return(tree)
  }

##############################################################################
#' localKappa
#'
#' A function that transforms branch lengths according to Kappa.
#' @param tree A tree of class phylo.
#' @param node A node number describing the clade to be transformed.
#' @param kappa The multiplier that rate is increased by in the clade specified by node.

# TODO(hfg): Add rescale option.

localKappa <- function(tree, node, kappa, rescale = TRUE) {
  descs <- getDescs(tree, node)
  descs <- c(node, descs)
  trans.edges <- which(tree$edge[ ,2] %in% descs)
  originalsum <- sum(tree$edge.length[trans.edges])
  res <- tree
  res$edge.length[trans.edges] <- tree$edge.length[trans.edges] ^ kappa
  
  if (rescale) {
    newsum <- sum(res$edge.length[trans.edges])
    ratio <- originalsum / newsum
    res$edge.length[trans.edges] <- res$edge.length[trans.edges] * ratio
  }
  
  return(res)
}

##############################################################################
#' localDelta
#'
#' A function that transforms a tree according to lambda, or a portion of a tree according to lambda.
#' @param tree A tree of class phylo.
#' @param node A node number describing the clade(s) to be transformed.
#' @param delta The value or values of lambda by which to transform the specified clade(s).
#' @name localDelta

# TODO(hfg): Add rescale option. Find the original root to tip length of the clade, divide it by the new root to tip length, and then multiply the new branch lengths by that.

localDelta <- function(tree, node, delta, rescale = TRUE) {
  descs <- getDescs(tree, node)
  tips <- tree$tip.label[descs[descs <= length(tree$tip.label)]]

  # Make a subtree by using drop.tip, and keep a copy of the edge matrix 
  # for comparison later on...

  subtree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% tips])

  n <- Ntip(subtree)
  hts <- data.frame(phytools::nodeHeights(subtree))
  colnames(hts) <- c("start", "end")
  T <- max(hts[,1])

  # Into heights I need to get the path length (the length from the root of the tree to that
  # node) and the branch length (which ought to be the end - the start)

  hts$t <- T - hts$end
  hts$bl <- hts$end - hts$start

  bls <- (hts$start + hts$bl) ^ delta - hts$start ^ delta

  subtree$edge.length <- bls  

  if (rescale) {
    scale <- T ^ delta
    subtree$edge.length <- (subtree$edge.length / scale) * T
  }

  tree$edge.length[tree$edge[ , 2] %in% descs] <- subtree$edge.length
  res <- tree
  return(res)
}

##############################################################################
#' localEB
#'
#' A function to transform a tree according to the exponential change/Early Burst model of Harmon et al 2010
#' @param tree A tree object of class phylo
#' @param node The node that the transformation should be applied at (includes all branches and nodes downstream of here).
#' @param a The EB parameter for the transformation.
#' @param rescale Rescale the tree after transformation? Defailts to TRUE.

localEB <- function(tree, node, a, rescale = TRUE) {
  print("Warning! This might be wrong (24/02/2017)")
  descs <- getDescs(tree, node)
  descs <- c(descs, node)
  trans.edges <- which(tree$edge[ ,1] %in% descs & tree$edge[ ,1] > length(tree$tip.label))
  bls <- tree$edge.length[trans.edges]
  times <- branching.times(tree)
  times <- times[names(times) %in% descs]
  maxt <- max(times)
  originalsum <- sum(tree$edge.length[trans.edges])
  res <- tree
  for (i in trans.edges) {
    branch <- res$edge.length[i]
    age <- times[which(names(times) == res$edge[i, 1])]
    t1 <- max(times) - age
    t2 <- t1 + branch
    res$edge.length[i] <- (exp(a * t2) - exp(a * t1))/(a)
  }
  
  if (rescale) {
    newsum <- sum(res$edge.length[trans.edges])
    ratio <- originalsum / newsum
    res$edge.length[trans.edges] <- res$edge.length[trans.edges] * ratio
  }
  
  return(res)
}

##############################################################################
#' treeTrans
#'
#' A function that transforms a tree, or a local part of a tree, according to lambda, kappa, delta, or a local rate. 
#' @param tree A tree of class phylo
#' @param param The transformation applied to the tree, either "lambda", "kappa", "delta", "rate" or "EB".
#' @param nodes A node number or vector of nodes describing the clade(s) to be transformed.
#' @param tips A vector, or list of vectors, of tip labels definining clade(s) to be transformed
#' @param value The value or vector of values to apply to the tree or parts of the tree. Order corresponds to the order of the elements of nodes or tips.
#' @param rescale Whether or not to rescale the tree to have the same root-to-tip length after transformation. Defaults to TRUE.
#' @export
#' @keywords tree transformation kappa lambda delta rates local rates local transformation
#' @examples
#' transTree(tree, param = "lambda", nodes = 52, value = 0)
#' transTree(tree, param = "rate", nodes = c(52, 91), value = 3)
#' transTree(tree, param = "delta", tips = list(c("dog", "cat", "moose"), c("frog", "salamander", "newt")), value = c(0.3, 2))

treeTrans <- function(tree, param, nodes = "root", tips = NULL, value, rescale = FALSE) {
  
  if (nodes == "root") {
    nodes <- length(tree$tip.label) + 1
  }

  if (is.null(nodes) & is.null(tips)) {
    stop("Must specify either node(s) or tips")
  }
  
  if (is.null(tips)) {
    
    if (length(nodes) != length(value)) {
      stop("Length of nodes and parameter values do not match.")
    }
    
    if (param == "lambda") {

      for (i in 1:length(nodes)) {
        tree <- localLambda(tree, nodes[i], value[i])
      }

    } else if (param == "kappa") {

      for (i in 1:length(nodes)) {
        tree <- localKappa(tree, nodes[i], value[i], rescale = rescale)
      }

    } else if (param == "delta") {

      for (i in 1:length(nodes)) {
        tree <- localDelta(tree, nodes[i], value[i], rescale = rescale)
      }

    } else if (param == "rate") {

      for (i in 1:length(nodes)) {
        tree <- localRate(tree, nodes[i], value[i])
      }
      
    } else if (param == "EB") {
    
      for (i in 1:length(nodes)) {
        tree <- localEB(tree, nodes[i], value[i], rescale = rescale)
      }
    
    }
  
  } else if (is.null(nodes)) {
    if (param == "lambda") {

      for (i in 1:length(tips)) {
        node <- getMRCA(tree, tips[[i]])
        tree <- localLambda(tree, node, value[i])
      }

    } else if (param == "kappa") {

      for (i in 1:length(tips)) {
        node <- getMRCA(tree, tips[[i]])
        tree <- localKappa(tree, node, value[i], rescale = rescale)
      }
    
    } else if (param == "delta") {

      for (i in 1:length(tips)) {
        node <- getMRCA(tree, tips[[i]])
        tree <- localDelta(tree, node, value[i], rescale = rescale)
      }
  
    } else if (param == "rate") {

      for (i in 1:length(tips)) {
        node <- getMRCA(tree, tips[[i]])
        tree <- localRate(tree, node, value[i])
      }

    }  else if (param == "EB") {
    
      for (i in 1:length(tips)) {
        node <- getMRCA(tree, tips[[i]])
        tree <- localEB(tree, nodes, value[i], rescale = rescale)
      }
    
    }
  
  
  }
  return(tree)
}


