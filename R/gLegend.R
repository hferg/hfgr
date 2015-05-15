#' gLegend
#'
#' A function to strip the legend from a ggplot object so that it can then be resued. For
#' example, you might want to plot 6 figures on one page with the same legend. This function
#' can take the legend from one of the plots, and then plot it later on after the 6 plots
#' are plotted with no legend. Normally used as part of ggplotSharedLegend.
#' HFG: From https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
#' @param a.gplot A ggplot object.
#' @export
#' @examples
#' leg <- gLegend(my_plot)
#' @import ggplot2

gLegend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
  }
  
