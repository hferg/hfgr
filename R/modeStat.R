#' modeStat
#' Calculate the mode.
#' @export

modeStat <- function(x) {
  z <- unique(x)
  x[which.max(tabulate(match(x, z)))]
}
