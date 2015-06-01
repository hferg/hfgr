#' massLh
#'
#' A function that will get the likelihoods for all the starter and true log files
#' in a directory, and return a data frame of those values for future comparison.
#' @param prefix The prefix of the log files. Needs to be followed, in the logfiles file names, with a number, then ".txt.starter.log" or ".txt.true.log"

massLh <- function(prefix) {
  files <- list.files()[grep(prefix, list.files())]
  starters <- files[grep("starter", files)]
  trues <- files[grep("true", files)]
  names <- gsub("\\..*", "", starters)
  res <- matrix(ncol = 4, nrow = length(starters))
  res[ ,1] <- names
  
  for (i in 1:nrow(res)) {
    raw.start <- readLines(starters[i])
    raw.true <- readLines(trues[i])
    start.output <- do.call(rbind, strsplit(raw.start[grep("\\bTree No\\b", raw.start):(length(raw.start) -1)], "\t"))
    true.output <- do.call(rbind, strsplit(raw.true[grep("\\bTree No\\b", raw.true):(length(raw.true) -1)], "\t"))
    res[i, 2] <- start.output[2, 2]
    res[i, 3] <- true.output[2, 2]
  }
  colnames(res) <- c("names", "startertreeLh", "truetreeLh", "diff")  
  res <- data.frame(res, stringsAsFactors = FALSE)
  res$startertreeLh <- as.numeric(res$startertreeLh)
  res$truetreeLh <- as.numeric(res$truetreeLh)
  res$diff <- res$truetreeLh - res$startertreeLh
  return(res)
}

