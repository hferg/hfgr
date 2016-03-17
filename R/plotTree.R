#' plotTree
#'
#' A silly little function that makes plotting a tree without tip labels and with
#' an axis easier and quicker.
#' @param tree An object of class phylo
#' @param tips Logical - show tip labels or not?
#' @param nodes TRUE shows all node labels, FALSE supresses node labels (default), or a vecotr of which nodes to highlight.
#' @param scale Logical - show the scale bar, or not?
#' @param cex Scaling factor for tip labels and node labels - universal.
#' @name plotTree
#' @export

plotTree <- function(tree, tips = FALSE, nodes = NULL, scale = TRUE, cex = 1) {
  plot(tree, show.tip.label = tips, cex = cex)
  if (!is.null(nodes)) {
    if (is.numeric(nodes)) {
      nodelabels(node = nodes, cex = cex)
    } else if (nodes == FALSE) {}
    else {
    nodelabels(cex = cex)
    }
  }
  if (scale) {
    axisPhylo()
  }
}
