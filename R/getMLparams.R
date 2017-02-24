#' getMLparams
#'
#' A function that will return the ML estimates from the logfile of an ML bayestraits run.
#' @param logfile The name of the logfile that needs to be examined, or vector of logfile names.
#' @name getMLparams
#' @export

getMLparams <- function(logfile) {
  if (length(logfile) == 1) {
    raw <- readLines(logfile)
    output <- do.call(rbind, strsplit(raw[grep("\\bTree No\\b", raw):length(raw)], "\t"))
    colnames(output) <- output[1, ]
    res <- output[c(2:nrow(output)), ]
  } else {
    res <- matrix(nrow = length(logfile), ncol = 5)
    for (i in 1:length(logfile)) {
      raw <- readLines(logfile[i])
      output <- do.call(rbind, strsplit(raw[grep("\\bTree No\\b", raw):length(raw)], "\t"))
      colnames(res) <- c("log", output[1, ])
      res[i, 1] <- logfile[i]
      res[i, c(2:ncol(res))] <- output[c(2:nrow(output)), ]
    }
  res <- data.frame(res)
  }
  return(res)
}

