#' loadRJ
#'
#' Returns the full mcmc object from a BayesTraits log file. This
#' is used inside plot functions and so on, but might be useful for
#' other MCMC manipulations and so on.
#' @param logfile The name of the logfile of the BayesTraits analysis.
#' @return A list containing the taxa translation table, all possible subtrees a scalar can occur on, and a data frame of the rj model configuration.
#' @export

loadRJ <- function(logfile) {

  raw <- readLines(logfile)
  rawhead <- strsplit(raw[1:(grep("\\bIt\\b", raw) -1)], "\t")
  rawtail <- strsplit(raw[grep("\\bIt\\b", raw):length(raw)], "\t")
  nms1 <- rawtail[[1]][1:7]
  nms2 <- rawtail[[1]][8:length(rawtail[[1]])]

  for (i in 1:length(rawtail)) {
    if (length(rawtail[[i]]) == 7) {
      names(rawtail[[i]]) <- nms1
    } else {
      len <- length(rawtail[[i]][8:length(rawtail[[i]])])
      end <- vector()

      for (j in 1:(len/4)) {
        end <- c(end, paste(nms2, j, sep = "_"))
      }
      nms <- c(nms1, end)
      names(rawtail[[i]]) <- nms
    }
  }
  
  # That has sorted out the RJ portion of the output.
  output <- do.call(smartBind, rawtail)
  output <- data.frame(output[2:nrow(output), ], stringsAsFactors = FALSE)

  # Now sort out the bit that describes the nodes, and defines the taxa.
  
  tipnum <- rawhead[[1]]
  taxatrans <- do.call(rbind, rawhead[c(1:tipnum+1)])
  
  subtreestart <- nrow(taxatrans) + 3
  subtrees <- rawhead[subtreestart:length(rawhead)]
  
  # Now name each element for each subtree.
  for (i in 1:length(subtrees)) {
    names(subtrees[[i]]) <- c(1:length(subtrees[[i]]))
  }
  
  subtrees <- do.call(smartBind, subtrees)
  subtrees <- data.frame(subtrees, stringsAsFactors = FALSE)
  colnames(subtrees)[c(1:2)] <- c("node", "bl")

  res <- list(taxatrans, subtrees, output)
  names(res) <- c("taxa", "subtrees", "rj_output")
  return(res)
}

