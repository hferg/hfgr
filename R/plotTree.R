#' plotTree
#'
#' A silly little function that makes plotting a tree without tip labels and with
#' an axis easier and quicker.
#' @param tree An object of class phylo
#' @param tips Logical - show tip labels or not?
#' @param nodes Logical - show nodel labels, or nor?
#' @param scale Logical - show the scale bar, or not?
#' @param cex Scaling factor for tip labels and node labels - universal.
#' @name plotTree
#' @export

plotTree <- function(tree, tips = FALSE, nodes = FALSE, scale = TRUE, cex = 1) {
  plot(tree, show.tip.label = tips, cex = cex)
  if (nodes) {
    nodelabels(cex = cex)
  }
  if (scale) {
    axisPhylo()
  }
}
