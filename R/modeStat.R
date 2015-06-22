#' modeStat
#' Calculate the mode.
#' @export

modeStat <- function(x, log = FALSE) {
  
  if (log == TRUE) {
    dens <- density(log(x))
    mode.tmp <- dens$x[which(dens$y == max(dens$y))]
    mode <- exp(mode.tmp)
  } else {
    dens <- density(x)
    mode <- dens$x[which(dens$y == max(dens$y))]
  }
  return(mode)
}

