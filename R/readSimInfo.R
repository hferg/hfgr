#' readSimInfo
#' A function that will read the simulation details from Andrews conDatSim function.
#' @param filename The filename of the simulation details file.
#' @param tree The tree the data were simulated over
#' @name readSimInfo
#' @export

readSimInfo <- function(filename, tree) {
  tree <- ladderize(tree)
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
    hds <- c(strsplit(cur[2], "\t")[[1]][c(1:8)], "ancNode", "descNode", "taxa")
    dat <- strsplit(cur[c(3:length(cur))], "\t")
    .res <- matrix(ncol = length(dat[[1]]) + 2, nrow = length(dat))  

    for (j in 1:length(dat)) {
      taxa <- strsplit(dat[[j]][9], ",")[[1]] 
      .res[j, c(1:8)] <- dat[[j]][c(1:8)]
      .res[j, c(9, 10)] <- tree$edge[which(tree$edge[ , 2] == mrca), ]
      .res[j, 11] <- dat[[j]][9]      
    }

    if (length(taxa) == 1) {
      mrca <- which(tree$tip.label == taxa)
    } else {
      mrca <- getMRCA(tree, taxa)
    }



    colnames(.res) <- hds
    res[[i]] <- .res
  }
  names(res) <- nms
  return(res)
}
