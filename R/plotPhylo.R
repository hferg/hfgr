#' plotPhylo
#'
#' A silly little function that makes plotting a tree without tip labels and with
#' an axis easier and quicker.
#' @param tree An object of class phylo
#' @param tips Logical - show tip labels or not?
#' @param nodes TRUE shows all node labels, FALSE supresses node labels (default), or a vecotr of which nodes to highlight.
#' @param scale Logical - show the scale bar, or not?
#' @param cex Scaling factor for tip labels and node labels - universal.
#' @param main A string that will be the title on the plot (if desired)
#' @param edge.cols A vector of edge colours, if you want the edges to be coloured.
#' @name plotPhylo
#' @export

plotPhylo <- function(tree, tips = FALSE, nodes = NULL, scale = TRUE, cex = 1, main = "",
  edge.cols = NULL) {
  plot(tree, show.tip.label = tips, cex = cex, main = main, edge.cols = edge.cols)
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
