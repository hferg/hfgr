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
  rawhead <- strsplit(raw[1:(grep("\\bIt*\\b", raw) -1)], "\t")
  rawtail <- strsplit(raw[grep("\\bIt*\\b", raw):length(raw)], "\t")
  nms1 <- rawtail[[1]][1:7]
  nms2 <- rawtail[[1]][8:length(rawtail[[1]])]
  nms2 <- gsub(" ", "", nms2)
  nms2 <- gsub("/", "", nms2)  

  for (i in 1:length(rawtail)) {
    if (length(rawtail[[i]]) == 7) {
      names(rawtail[[i]]) <- nms1
    } else {
      len <- length(rawtail[[i]][8:length(rawtail[[i]])])
      end <- vector(mode = "character", length = len)

      for (j in 1:length(end)) {
        end[j] <- paste(nms2, j, sep = "_")
      }
      nms <- c(nms1, end)
      names(rawtail[[i]]) <- nms
    }
  }

  tipnum <- rawhead[[1]]
  taxatrans <- do.call(rbind, rawhead[c(1:tipnum+1)])  
  subtreestart <- nrow(taxatrans) + 3
  subtrees <- rawhead[subtreestart:length(rawhead)]
  
  for (i in 1:length(subtrees)) {
    names(subtrees[[i]]) <- c(1:length(subtrees[[i]]))
  }
  
  output <- do.call(smartBind, rawtail)
  output <- data.frame(output[2:nrow(output), ], stringsAsFactors = FALSE)
  subtrees <- do.call(smartBind, subtrees)
  subtrees <- data.frame(subtrees, stringsAsFactors = FALSE)
  colnames(subtrees)[c(1:2)] <- c("node", "bl")

  res <- list(taxatrans, subtrees, output)
  names(res) <- c("taxa", "subtrees", "rj_output")
  return(res)
}

