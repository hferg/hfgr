# Here is a function for you. pp is the output of the post-processor, tree is the timetree, and threshold is the decimal of the threshold you want.
treez4chris <- function(pp, threshold, measure = "mean") {
  if (measure == "mean") {
    ms <- "meanBL"
  } else if (measure == "median") {
    ms <- "medianBL"
  }

  tree <- pp$meantree
  tree$edge.length <- pp$data$orgBL[-1]
  perscl <- apply(pp$scalars$rates, 1, function(x) sum(x != 1)) / pp$niter
  scaled_branches <- as.numeric(names(which(perscl[-1] > threshold)))
  tree$edge.length[tree$edge[ , 2] %in% scaled_branches] <- pp$data[pp$data$descNode %in% scaled_branches, ms]
  return(tree)
}
