#'
#' cladeChangeSim
#' Simulates the expected change along all branches of a defined clade.
#' @name cladeChangeSim
#' @export

cladeChangeSim <- function(tree, reftree = NULL, node, mode = "BM", directional = "non", startval = 0, sigsq, 
  plot = FALSE, increment = NULL, param.value = NULL, segments = FALSE) {
    # Identify all edges - getDesTips_edge returns tip LABELS descendent from the branch (edge)
    # As such, this gives you all of the branches in the clade.
    # HFG CHANGE HERE
    descs <- getDescs(tree, node)
    tips <- descs[which(descs <= length(tree$tip.label))]
    internal <- descs[which(descs > length(tree$tip.label))]
   
 
    allbranches <- c(which((tree$edge[ , 2] == node)), which.edge(tree, tree$tip.label[tips]))
    
    res <- vector(mode = "list", length = length(allbranches))
    names(res) <- allbranches

    if (is.null(increment)) {
      increment <- max(tree$edge.length[allbranches]) / 500
    }

    if (is.null(reftree)) {
      for (i in allbranches) {
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
        if (segments) {
          sim[ , 1] <- sim[ , 1]+ getPL(tree, node = tree$edge[i, 1])
          .tmp <- matrix(ncol = 5, nrow = (nrow(sim) - 1))
          colnames(.tmp) <- c("start", "end", "ancstate", "trait", "branch")
          .tmp[ , "start"] <- sim[c(1:nrow(sim)) - 1, 1]
          .tmp[ , "end"] <- sim[c(2:nrow(sim)) , 1] 
          .tmp[ , "ancstate"] <- sim[c(1:nrow(sim)) - 1, 2]
          .tmp[ , "trait"] <- sim[c(2:nrow(sim)), 2]
          .tmp[ , "branch"] <- i
          sim <- .tmp
          res[[which(names(res) == i)]] <- sim
        } else {
          colnames(sim) <- c("time", "trait")
          res[[which(names(res) == i)]] <- sim
          res[[which(names(res) == i)]][ , 1] <- res[[which(names(res) == i)]][ , 1] + getPL(tree, node = tree$edge[i, 1])          
        }
      }
    } else {
      for (i in allbranches) {
        bl <- tree$edge.length[i]
        refbl <- reftree$edge.length[i]
        blrat <- refbl / bl
        leading <- which(tree$edge[ , 2] == tree$edge[i, 1])

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
        if (segments) {
          sim[ , 1] <- sim[ , 1]+ getPL(tree, node = tree$edge[i, 1])
          .tmp <- matrix(ncol = 5, nrow = (nrow(sim) - 1))
          colnames(.tmp) <- c("start", "end", "ancstate", "trait", "branch")
          .tmp[ , "start"] <- sim[c(1:nrow(sim)) - 1, 1]
          .tmp[ , "end"] <- sim[c(2:nrow(sim)) , 1]
          .tmp[ , "ancstate"] <- sim[c(1:nrow(sim)) - 1, 2]
          .tmp[ , "trait"] <- sim[c(2:nrow(sim)), 2]
          .tmp[ , "branch"] <- i
          sim <- .tmp
          res[[which(names(res) == i)]] <- sim
        } else {
          colnames(sim) <- c("time", "trait")
          res[[which(names(res) == i)]] <- sim
          res[[which(names(res) == i)]][ , 1] <- res[[which(names(res) == i)]][ , 1] + getPL(tree, node = tree$edge[i, 1])          
        }
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
    test <- melt(res)
    part1 <- test[test$Var2 == "time", ]
    part2 <- test[test$Var2 == "trait", ]
    res <- data.frame(time = part1$value, trait = part2$value, branch = part1$L1)    
    if (plot) {
      p <- ggplot(res, aes(x = time, y = trait, group = branch))
      print(p + geom_line())
    }
  }
  return(res)
}
