#' ggplotShareLeg
#'
#' A function to plot multiple ggplot objects with a shared legend.
#' HFG: from https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
#' @param ... Some number of ggplot objects.
#' @return A multiple panel plot of plots with a common legend.
#' @export
#' @examples
#' dsamp <- diamonds[sample(nrow(diamonds), 1000), ]
#' p1 <- qplot(carat, price, data=dsamp, colour=clarity)
#' p2 <- qplot(cut, price, data=dsamp, colour=clarity)
#' p3 <- qplot(color, price, data=dsamp, colour=clarity)
#' p4 <- qplot(depth, price, data=dsamp, colour=clarity)
#' ggplotShareLeg(p1, p2, p3, p4)
#' @import ggplot2
#' @import gridExtra

ggplotShareLeg <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
        x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

