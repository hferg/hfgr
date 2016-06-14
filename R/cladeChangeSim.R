#'
#' cladeChangeSim
#' Simulates the expected change along all branches of a defined clade.
#' @param tree The time tree over which to simulate change.
#' @param reftree A reference tree that may have divergent branch lengths compared the the time tree (for example, a rate scaled tree)
#' @param node The node from which the simulation starts (defaults to root)
#' @param plot Plot the sim to screen after running?
#' @param lead Include the leading branch in the simulation?
#' @name cladeChangeSim
#' @export

cladeChangeSim <- function(tree, reftree = NULL, node = "root", mode = "BM", directional = "non", 
    startval = 0, sigsq, plot = FALSE, increment = NULL, segments = FALSE, lead = FALSE) {
    # Identify all edges - getDesTips_edge returns tip LABELS descendent from the branch (edge)
    # As such, this gives you all of the branches in the clade.
    # HFG CHANGE HERE
    
  if (node == "root") {
    node <- length(tree$tip.label) + 1
  }

  descs <- getDescs(tree, node)
  tips <- descs[which(descs <= length(tree$tip.label))]
  internal <- descs[which(descs > length(tree$tip.label))]

  allbranches <- findBranches(tree, node, tail = lead)
  
  res <- vector(mode = "list", length = length(allbranches))
  names(res) <- allbranches

  if (is.null(increment)) {
    increment <- max(tree$edge.length[allbranches]) / 500
  }

  for (i in allbranches) {
    if (is.null(reftree)) {
      bl <- tree$edge.length[i]
      blrat <- 1
      leading <- which(tree$edge[ , 2] == tree$edge[i, 1])
    } else {
      bl <- tree$edge.length[i]
      refbl <- reftree$edge.length[i]
      blrat <- refbl / bl
      leading <- which(tree$edge[ , 2] == tree$edge[i, 1])
    }

    if (length(leading) == 0) {
      strt <- startval
    } else if (!leading %in% allbranches) {
      strt <- startval
    } else {
      tmp <- res[[which(names(res) == leading)]]
      strt <- tmp[nrow(tmp), "trait"]
    }
    
    sim <- branchChangeSim(branchtime = bl, direction = directional, startval = strt, var = (sigsq * blrat), 
        increment = increment)
    .sim <- matrix(ncol = 3, nrow = nrow(sim))
    .sim[ , 3] <- tree$edge[i, 2]
    .sim[ , c(1, 2)] <- sim[ , c(1, 2)]
    sim <- .sim
    if (segments) {
      sim[ , 1] <- sim[ , 1]+ getPL(tree, node = tree$edge[i, 1])
      .tmp <- matrix(ncol = 6, nrow = (nrow(sim) - 1))
      colnames(.tmp) <- c("start", "end", "ancstate", "trait", "branch", "descNode")
      .tmp[ , "start"] <- sim[c(1:nrow(sim)) - 1, 1]
      .tmp[ , "end"] <- sim[c(2:nrow(sim)) , 1]
      .tmp[ , "ancstate"] <- sim[c(1:nrow(sim)) - 1, 2]
      .tmp[ , "trait"] <- sim[c(2:nrow(sim)), 2]
      .tmp[ , "branch"] <- i
      .tmp[ , "descNode"] <- tree$edge[i , 2]
      sim <- .tmp
      res[[which(names(res) == i)]] <- sim
    } else {
      colnames(sim) <- c("time", "trait", "descNode")
      res[[which(names(res) == i)]] <- sim
      res[[which(names(res) == i)]][ , 1] <- res[[which(names(res) == i)]][ , 1] + getPL(tree, node = tree$edge[i, 1])          
    }
  }

  if (segments) {
    res <- do.call(rbind, res)
    res <- as.data.frame(res)
    if (plot) {
      p <- ggplot(res_seg)
      t <- geom_segment(data = res_seg, aes(x = start, xend = end, y = ancstate, yend = trait), alpha = 0.8)
      print(p + t)
    }    
  } else {
    res <- melt(res)
    time <- res[res$Var2 == "time", ]
    trait <- res[res$Var2 == "trait", ]
    descNode <- res[res$Var2 == "descNode", ]
    res <- data.frame(time = time$value, trait = trait$value, branch = time$L1, descNode = descNode$value)    
    if (plot) {
      p <- ggplot(res, aes(x = time, y = trait, group = branch))
      print(p + geom_line())
    }
  }
  return(res)
}
