#' treeTrans
#'
#' A function that transforms a tree, or a local part of a tree, according to lambda, kappa, delta, or a local rate. 
#' @param tree A tree of class phylo
#' @param param The transformation applied to the tree, either "lambda", "kappa", "delta", or "rate"
#' @param nodes A node number or vector of nodes describing the clade(s) to be transformed.
#' @param tips A vector, or list of vectors, of tip labels definining clade(s) to be transformed
#' @param value The value or vector of values to apply to the tree or parts of the tree. Order corresponds to the order of the elements of nodes or tips.
#' @export
#' @keywords tree transformation kappa lambda delta rates local rates local transformation
#' @examples
#' transTree(tree, param = "lambda", nodes = 52, value = 0)
#' transTree(tree, param = "rate", nodes = c(52, 91), value = 3)
#' transTree(tree, param = "delta", tips = list(c("dog", "cat", "moose"), c("frog", "salamander", "newt")), value = c(0.3, 2))

treeTrans <- function(tree, param, nodes = NULL, tips = NULL, value, rescale = TRUE) {
  
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

    }
  
  }
  return(tree)
}

