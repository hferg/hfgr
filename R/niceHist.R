#' niceHist
#' 
#' A function to plot a nice histogram - one that uses ggplot2
#' but has sensible break behaviour in R. Produces a single
#' histogram but is mostly called internally for things like
#' plotPosts.
#' @param parameter The a vecotr of values that you want a histogram of.
#' @return A histogram with nice ggplot aesthetics, and good bin widths.
#' @export

niceHist <- function(parameter, fill = "dodgerblue", breaks = "scott") {
  dat <- parameter[!is.na(parameter)]
  
  if (breaks == "scott") {
    bwidth <- 3.5 * sd(dat) * length(dat) ^ -(1/3)
  } else if (breaks == "fandd") {
    bwidth <- 2 * IQR(dat) * length(dat) ^ -(1/3)
  } else {
    bwidth = (max(dat) - min(dat))/30
  }
    
  if (fill == "count") {
    ggplot(data.frame(parameter = dat), aes(x = parameter, fill = ..count..)) +
      geom_histogram(color = "darkgray", binwidth = bwidth)
  } else {
    ggplot(data.frame(parameter = dat), aes(x = parameter)) +
      geom_histogram(color = "darkgray", fill = fill, binwidth = bwidth)
  }
}

