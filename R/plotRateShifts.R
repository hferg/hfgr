#' plotRateShifts
#' A simple function to scale the tree scaled by the mean, median or mode scalars of branches
#' that are scaled over a certain threshold.
#' @param pp An output object from rjpp
#' @param threshold The threshold over which branches are scaled (in decimal form)
#' @param measure "mean", "median" or "mode" - how to summarise the scalars that are used.
#' @name plotRateShifts
#' @export

plotRateShifts <- function(pp, threshold, measure = "mean") {
  if (measure == "mean") {
    ms <- "meanBL"
  } else if (measure == "median") {
    ms <- "medianBL"
  } else if (measure == "mode") {
    ms <- "modeBL"
  }

  tree <- pp$meantree
  tree$edge.length <- pp$data$orgBL[-1]
  perscl <- apply(pp$scalars$rates, 1, function(x) sum(x != 1)) / pp$niter
  scaled_branches <- as.numeric(names(which(perscl[-1] > threshold)))
  tree$edge.length[tree$edge[ , 2] %in% scaled_branches] <- pp$data[pp$data$descNode %in% scaled_branches, ms]
  return(tree)
}
