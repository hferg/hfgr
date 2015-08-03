#' loadRJdisc
#'
#' Returns the full mcmc object from a BayesTraits discrete model log file. This
#' is used inside plot functions and so on, but might be useful for
#' other MCMC manipulations and so on.
#' @param logfile The name of the logfile of the BayesTraits analysis.
#' @param burnin burning
#' @param thinning Thin the sample by ever n sample.
#' @return A list containing the taxa translation table, all possible subtrees a scalar can occur on, and a data frame of the rj model configuration.
#' @export
#' @name loadRJdisc

loadRJdisc <- function(logfile, burnin = 0, thinning = 1) {

  raw <- readLines(logfile)
  rawhead <- strsplit(raw[1:(grep("\\bIteration*\\b", raw) -1)], "\t")
  rawtail <- strsplit(raw[grep("\\bIteration*\\b", raw):length(raw)], "\t")
  nms1 <- rawtail[[1]][1:7]
  nms2 <- rawtail[[1]][8:length(rawtail[[1]])]
  nms2 <- gsub(" ", "", nms2)
  nms2 <- gsub("/", "", nms2)  

  rjout <- do.call(rbind, rawtail[c(2:length(rawtail))])
  rjout <- data.frame(rjout)
  colnames(rjout) <- c(nms1, nms2)

  # Convert all columns except for model string and
  cols <- colnames(rjout)
  for (i in colnames(rjout)) {
    if (i != "Harmonic Mean" & i != "Model string") {
      rjout[ , i] <- as.numeric(as.character(rjout[ , i]))
    }
  }
  rjout <- rjout[seq.int(burnin, nrow(rjout), thinning), ]
  return(rjout)
}
