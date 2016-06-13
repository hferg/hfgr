#' simOverTrees
#' This function will simulate the expected change over a vector of named trees. It calls cladeChangeSim a bunch of times.
#' 
#' @param tree The original tree.
#' @param reftrees A vector of trees that differ to the original tree on which to simulate the change.
#' @param nsim The number of simulations per tree.
#' @param sigsq A single sigma squared for the simulation, or a vector of them in the order of the trees.
#' @param node A node from which to conduct the simulation.
#' @param direction Positive selection ("pos"), negative selection ("neg"), or non directional ("non")?
#' @param increment The size of the time interval each branch is broken down into, if NULL is calculated automatically, but may not be comparable between runs.
#' @name simOverTrees
#' @export

simOverTrees <- function(tree, reftrees, nsim, node, sigsq, direction, increment = NULL) {
  
  if (node == "root") {
    node <- length(tree$tip.label) + 1
  } else {
    node = node
  }

  nsim <- nsim
  res <- vector(mode = "list", length = nsim)
  if (!is.null(names(reftrees))) {
    nms <- names(reftrees)
  } else {
    nms <- c(1:length(reftrees))
  }

  if (length(sigsq) == 1) {
    sigsq <- rep(sigsq, length(reftrees))
  }
  for (i in 1:nsim) {
    tmp_res <- vector(mode = "list", length = length(reftrees))
    names(tmp_res) <- nms
    for (j in 1:length(reftrees)) {
      tmp_res[[j]] <- cladeChangeSim(tree, reftree = reftrees[[j]], node = node, sigsq = sigsq[j], 
        directional = direction, increment = increment)
      tmp_res[[j]]$sim <- i
    }

    dat_tmp <- data.frame(tmp_res)
    dat <- dat_tmp[ , grep("trait", colnames(dat_tmp))]
    dat$time <- dat_tmp[ , grep("time", colnames(dat_tmp))[1]]
    dat$branch <- dat_tmp[ , grep("branch", colnames(dat_tmp))[1]]
    dat$sim <- dat_tmp[ , grep("sim", colnames(dat_tmp))[1]]
    res[[i]] <- dat
  }
  res <- do.call(rbind, res)
  return(res)
}
