#' testDelta
#' 
#' A function to test the output of Andrew's BM model to see if delta
#' has or has not been simulated using geiger to test the likelihood.
#' @param tree The tree the data was simulated on
#' @param node The node the delta is supposed to be at
#' @param datafile The name of the datafile the simulation outputs.
#' @export
#' @name testDelta

testDelta <- function(tree, datafile, node = "root") {
  if (node == "root") {
    node <- length(tree$tip.label) + 1
  }

  d <- read.table(datafile)
  tips <- getTipNames(tree, node)
  subtree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% tips])
  res <- matrix(ncol = 5, nrow = length(c(2:ncol(d))))
  colnames(res) <- c("bm", "bmss", "delta", "deltaval", "delss")
  
  for (i in 2:ncol(d)) {
    subdata <- d[ , i]
    names(subdata) <- d[ , 1]
    subdata <- subdata[names(subdata) %in% tips]
    bm <- fitContinuous(subtree, subdata, model = "BM")
    del <- fitContinuous(subtree, subdata, model = "delta")
    res[i-1, "bm"] <- bm$opt$lnL
    res[i-1, "bmss"] <- bm$opt$sigsq
    res[i-1, "delta"] <- del$opt$lnL
    res[i-1, "deltaval"] <- del$opt$delta
    res[i-1, "delss"] <- del$opt$sigsq
  }
  rownames(res) <- c(1:(ncol(d) - 1))
  return(res)
}
