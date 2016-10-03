#' posthocZ
#' 
#' Simulates a trait over the posterior sample of trees, and then calculates z scores.
#' @param tree The time tree
#' @param rjtrees The posterior of trees in multiphylo form, or the filename of the posterior.
#' @param nsim The number of simulations per tree in the posterior.
#' @export

posthocZ <- function(tree, rjtrees, nsim = 100) {
  if (class(rjtrees) == "multiPhylo") {
    posterior <- rjtrees
  } else {
    posterior <- read.nexus(rjtrees)
  }

  n <- rep(1:length(posterior), nsim)

  postsims <- do.call(rbind, lapply(1:length(n), function(x) rTraitCont(posterior[[n[x]]], sigsq = 1, model = "BM")))
  timesims <- do.call(rbind, lapply(1:length(n), function(x) rTraitCont(tree, sigsq = 1, model = "BM")))

  postsd <- apply(postsims, 2, sd)
  timesd <- apply(timesims, 2, sd)

  res <- data.frame(tip = names(postsd), sd = postsd, z = (postsd - mean(timesd)) / (sd(postsd) * sqrt((length(postsd) - 1) / length(postsd))))
  return(res)
}
