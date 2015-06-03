#' readPP
#'
#' A function that reads the output of the variable rates post-processor (VarPP-v1). This is the output that VarPP-v1 prints to screen, so that needs to be written to a text file first.
#' @param logfile The recorded output of VarPP-v1
#' @export

readPP <- function(logfile) {
  raw <- readLines(logfile)
  split <- strsplit(raw, "\t")
  nms <- split[[1]]
  nms <- gsub(" ", "_", nms)
  pp <- data.frame(do.call(rbind, split[c(2:length(split))]), stringsAsFactors = FALSE)
  colnames(pp) <- nms
  
  for (i in 1:8) {
    pp[ ,i] <- as.numeric(pp[ ,i])
  }
  
  return(pp)
}

