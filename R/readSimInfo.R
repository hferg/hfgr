#' readSimInfo
#' A function that will read the simulation details from Andrews conDatSim function.
#' @param filename The filename of the simulation details file.
#' @param tree The tree the data were simulated over
#' @name readSimInfo
#' @export

readSimInfo <- function(filename, tree) {
  raw <- readLines(filename)
  nsims <- length(grep("==", raw))
  simpos <- grep("==", raw)

  res <- vector(mode = "list", length = nsims)
  nms <- vector(mode = "character", length = nsims)

  for (i in 1:length(simpos)) {
    if (i != length(simpos)) {
      cur <- raw[c(simpos[i]:(simpos[i + 1] - 1))]
    } else {
      cur <- raw[c(simpos[i]:length(raw))]
    }

    nms[i] <- paste0(strsplit(gsub("=", "", cur[1]), "\t")[[1]], collapse = "")
    hds <- strsplit(cur[2], "\t")[[1]]
    dat <- strsplit(cur[c(3:length(cur))], "\t")
    .res <- matrix(ncol = length(dat[[1]]), nrow = length(dat))  

    for (j in 1:length(dat)) {
      .res[j, ] <- dat[[j]]
    }

    colnames(.res) <- hds
    res[[i]] <- .res
  }
  names(res) <- nms
  return(res)
}
