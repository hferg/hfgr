#'
#' cladeChangeSim
#' Simulates the expected change along all branches of a defined clade.
#' @name cladeChangeSim
#' @export

cladeChangeSim <- function(tree, node, directional = "non", startval = 0, sigsq, plot = FALSE) {
    # Identify all edges - getDesTips_edge returns tip LABELS descendent from the branch (edge)
    # As such, this gives you all of the branches in the clade.
    # HFG CHANGE HERE
    descs <- getDescs(tree, node)
    tips <- descs[which(descs <= length(tree$tip.label))]
    internal <- descs[which(descs > length(tree$tip.label))]
    allbranches <- c(which((tree$edge[ , 2] == node)), which.edge(tree, tree$tip.label[tips]))
    res <- vector(mode = "list", length = length(allbranches))
    names(res) <- allbranches

    for (i in allbranches) {
      increment <- max(tree$edge.length[allbranches]) / 500
      bl <- tree$edge.length[i]
      leading <- which(tree$edge[ , 2] == tree$edge[i, 1])

      if (length(leading) == 0) {
        strt <- startval
      } else if (!leading %in% allbranches) {
        strt <- startval
      } else {
        tmp <- res[[which(names(res) == leading)]]
        strt <- tmp[nrow(tmp), 2]
      }
      
      sim <- branchChangeSim(branchtime = bl, direction = directional, startval = strt, var = sigsq, increment = increment)      
      res[[which(names(res) == i)]] <- sim
      res[[which(names(res) == i)]][ , 1] <- res[[which(names(res) == i)]][ , 1] + getPL(tree, node = tree$edge[i, 1])
    }
  test <- melt(res)
  part1 <- test[test$Var2 == 1, ]
  part2 <- test[test$Var2 == 2, ]
  res <- data.frame(time = part1$value, trait = part2$value, branch = part1$L1)

  if (plot) {
    p <- ggplot(res, aes(x = time, y = trait, group = branch))
    print(p + geom_line())
  }

  return(res)
}
