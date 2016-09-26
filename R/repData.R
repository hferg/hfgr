#' repData
#' 
#' This function will write a vector of data points to the working directory n times, with a given prefix and the number of the replicate as the suffix. Ideally, the names of the vector are the tip labels of the tree!
#' @param dat The data file to be replicated.
#' @param n The number of times to write the file to the disk, or a vector of suffixes that the data is replicated using.
#' @param prefix The prefix for each of the copies of the data file.
#' @export

repData <- function(dat, n, prefix) {
  
  if (!is.vector(dat)) {
    stop("Requires a named vector of data.")
  }

  if (is.null(names(dat))) {
    stop("Requires a named vector of data.")
  }
  
  if (is.vector(n)) {
    x <- n
  } else {
    x <- 1:n
  }
  
  invisible(sapply(x, function(x) write.table(dat, file = paste0(prefix, "_", x, ".txt"), quote = FALSE, col.names = FALSE)))

}
