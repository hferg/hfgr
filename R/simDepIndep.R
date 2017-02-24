#' simDepIndep
#'
#' A function that simulates either dependent or independent evolution between two binary
#' character over a tree, and then takes some number of clades of a specific size, and simulates
#' dependent or independent evolution (whichever the rest of the tree isn't) within those clades.
#' The tree and data are simulated, and the Q matrices for the two traits are fixed for both 
#' and independent evolution (at the moment).
#' @export
#' @name simDepIndep
#' @param itts The number of trees and datasets to generate.
#' @param treesize The size, in terminal taxa, of the tree to be simulated.
#' @param mintax The minimum size for the clades that vary from the background model.
#' @param maxtax The maximum size for the clades that vary from the background model.
#' @param clades The number of clades to simulate the opposite model in.
#' @param base_mode The model for the whole tree that the clades vary from. E.g. if "dependent" then the tree will be a dependent model, with n clades evolved as independent. If "independent" than vice versa.
#' @param base_rate The model defaults to a "base transition rate" of 2 changes per mean path length, this changes that.

simDepIndep <- function(itts, treesize, mintax, maxtax, clades, base_mode, base_rate = 2) {
  trees <- vector(mode = "list", length = itts)
  for (i in 1:itts) {
    candidates <- 0
    while (length(candidates) <= clades) {
      tree <- sim.bd.taxa(treesize, numbsim = 1, lambda = 0.2, mu = 0.01, complete = FALSE)[[1]]
      tree <- ladderize(tree)
      tree$edge.length <- tree$edge.length / max(nodeHeights(tree))
      nodes <- matrix(nrow = nrow(tree$edge), ncol = 2)
      colnames(nodes) <- c("Node", "nTips")
      nodes[ , 1] <- tree$edge[ , 2]

      meanpl <- mean(nodeHeights(tree)[ , 2])
      # btrans is the base transition rate
      btrans <- round(base_rate / meanpl, 1)

      for (k in 1:nrow(nodes)) {
        nodes[k, 2] <- sum(getDescs(tree, nodes[k, 1]) < length(tree$tip.label)) 
      }

      # Find two clades that have ~ 30 tips.
      candidates <- nodes[nodes[ , "nTips"] >= mintax & nodes[ , "nTips"] <= maxtax, "Node"]

      ## I need to make a way to automate finding the clades of interst, and ind oing so avoid the
      # situation where it finds two nodes that are nested within each other. The first thing is to
      # find nodes from the list of candidates that are not nested, and then pick at random from them.

      change_nodes <- NULL
      pairs <- list()

      # If I just find the nested pairs, then pick one from each, then discard them from the candidates,
      # I should be alright...
      for (k in candidates) {
        descs_focus <- getDescs(tree, k)
        remainder <- candidates[candidates != k]

        # If there are nested nodes, put them together into a different object
        # in order to later pick a single one of them at random.
        for (j in remainder) {
          descs_remainder <- getDescs(tree, j)
          if (any(descs_focus %in% descs_remainder)) {
            pairs[[length(pairs) + 1]] <- c(k, j)
          }
        }

      }

      if (length(pairs) > 0) {

        for (k in 1:length(pairs)) {
          pairs[[k]] <- sort(pairs[[k]])
        }

        pairs <- matrix(unlist(unique(pairs)),,2, byrow = TRUE)
        odds <- unique(as.vector(pairs) [sapply(as.vector(pairs), function(x) sum(pairs == x) >1)])
        pairs_n_odds <- c(pairs[!c(pairs) %in% odds], odds)

        pairs_selection <- NULL
        odds_selection <- NULL

        if (length(odds) > 0) {
          odds_selection <- sample(odds, 1)
        }

        for (k in 1:nrow(pairs)) {
          select <- sample(pairs[k, ], 1)
          if (!select %in% odds)  
            pairs_selection <- c(pairs_selection, select)
        }

        candidates <- c(candidates[!candidates %in% pairs_n_odds], odds_selection, pairs_selection)
      }
    }

    change_nodes <- sample(candidates, clades)
    change_taxa <-  vector(mode = "list", length = length(change_nodes))
    for (k in 1:length(change_taxa)) {
      tips <- getDescs(tree, change_nodes[[k]])
      tips <- tips[tips <= length(tree$tip.label)]
      tips <- tree$tip.label[tips]
      change_taxa[[k]] <- tips
    }

    # Now split off the subtrees for simulation.
    changetrees <- vector(mode = "list", length = clades)
    for (k in 1:clades) {
      tips <- getDescs(tree, change_nodes[k])
      tips <- tips[tips  < length(tree$tip.label)]
      tips <- tree$tip.label[tips]
      changetrees[[k]] <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% tips])
    }

    # now simulate the data.
    if (base_mode == "dependent") {
      Q <- matrix(c(0,          btrans / 2, btrans / 2, 0,
                    btrans * 2, 0,          0,          btrans * 2,
                    btrans * 2, 0,          0,          btrans * 2,
                    0,          btrans / 2, btrans / 2, 0), 
          4, 4, byrow=TRUE)
      # Name it and set the diagonal
      rownames(Q) <- colnames(Q) <- c("aa", "ab", "ba", "bb")
      diag(Q) <- -rowSums(Q)
      # Now simulate the history for the tree (this is basically a 4 state single trait)
      tt <- sim.history(tree, Q)
      # And now split the 4 state single trait into two 2-state traits.
      t1 <- mergeMappedStates(tt, c("aa", "ab"), "a")
      t1 <- mergeMappedStates(t1, c("ba", "bb"), "b")
      t2 <- mergeMappedStates(tt, c("aa", "ba"), "a")
      t2 <- mergeMappedStates(t2, c("ab", "bb"), "b")
      t1$states <- getStates(t1, "tips")
      t2$states <- getStates(t2, "tips")

      # Now turn those data into a table for use in the model.
      dat <- as.data.frame(cbind(t1$states, t2$states))

      # Now find the root state for trait 1 for each of the change_nodes
      for (k in 1:length(changetrees)) {
        Q_ind <- matrix(c(-btrans, btrans, btrans, -btrans), 2)
        rownames(Q_ind) <- colnames(Q_ind) <- letters[1:2]
        # TRAIT 1      
        root <- names(t1$maps[[which(t1$edge[ , 2] == change_nodes[k])]])
        root <- root[length(root)]

        tmp <- sim.history(changetrees[[k]], Q_ind, anc = root)
        tmp <- tmp$states
        dat[rownames(dat) %in% names(tmp), "V1"] <- tmp

        # TRAIT 2
        root <- names(t2$maps[[which(t2$edge[ , 2] == change_nodes[k])]])
        root <- root[length(root)]

        tmp <- sim.history(changetrees[[k]], Q_ind, anc = root)
        tmp <- tmp$states

        dat[rownames(dat) %in% names(tmp), "V2"] <- tmp
      }
    
    } else if (base_mode == "independent") {

        Q_ind <- matrix(c(-btrans, btrans, btrans, -btrans), 2)
        rownames(Q_ind) <- colnames(Q_ind) <- letters[1:2]
        t1 <- sim.history(tree, Q_ind)
        t2 <- sim.history(tree, Q_ind)
        dat <- as.data.frame(cbind(t1$states, t2$states))

        for (k in 1:length(changetrees)) {
        
          Q_dep <- matrix(c(0,          btrans / 2, btrans / 2, 0,
                    btrans * 2, 0,          0,          btrans * 2,
                    btrans * 2, 0,          0,          btrans * 2,
                    0,          btrans / 2, btrans / 2, 0), 
          4, 4, byrow=TRUE)
          # Name it and set the diagonal
          rownames(Q_dep) <- colnames(Q_dep) <- c("aa", "ab", "ba", "bb")
          diag(Q_dep) <- -rowSums(Q_dep)
          # Since this is dependent the two traits are estimated simultaneously.
          # Establish root state for the subtree. This needs to be for both characters,
          # i.e. aa or bb, or something.

          root1 <- names(t1$maps[[which(t1$edge[ , 2] == change_nodes[k])]])
          root2 <- names(t2$maps[[which(t2$edge[ , 2] == change_nodes[k])]])
          root <- c(root1[length(root1)], root2[length(root2)])
          root <- paste0(root[1], root[2])

          # Now simulate some dependent data for the current changetree based on this
          # root node.

          tts <- sim.history(changetrees[[k]], Q_dep, anc = root)
          
          if (any(tts$states == "aa")) {
            t11 <- mergeMappedStates(tts, "aa", "a")
          }

          if (any(t11$states == "ab")) {
            t11 <- mergeMappedStates(t11, "ab", "a")
          }
          
          if (any(t11$states == "ba")) {
            t11 <- mergeMappedStates(t11, "ba", "b")
          }

          if (any(t11$states == "bb")) {
            t11 <- mergeMappedStates(t11, "bb", "b")
          }
          
          if (any(tts$states == "aa")) {
            t22 <- mergeMappedStates(tts, "aa", "a")
          }

          if (any(t22$states == "ba")) {
            t22 <- mergeMappedStates(t22, "ba", "a")
          }
          
          if (any(t22$states == "ab")) {
            t22 <- mergeMappedStates(t22, "ab", "b")
          }

          if (any(t22$states == "bb")) {
            t22 <- mergeMappedStates(t22, "bb", "b")
          }
          

          
          t11$states <- getStates(t11, "tips")
          t22$states <- getStates(t22, "tips")

          tmp <- as.data.frame(cbind(t11$states, t22$states))
          # Now merge these states back into the original data.
          dat[rownames(dat) %in% rownames(tmp), "V1"] <- tmp$V1
          dat[rownames(dat) %in% rownames(tmp), "V2"] <- tmp$V2
        }

    }

    # Now make dat numeric - swap a for 0 and b for 1.
    dat$V1 <- gsub("a", 0, dat$V1)
    dat$V2 <- gsub("a", 0, dat$V2)
    dat$V1 <- gsub("b", 1, dat$V1)
    dat$V2 <- gsub("b", 1, dat$V2)
    dat$V1 <- as.numeric(dat$V1)
    dat$V2 <- as.numeric(dat$V2)
    colnames(dat) <- c("trait1", "trait2")

    trees[[i]] <- list(tree = tree, data = dat, changed_taxa = change_taxa)
  }
    
  return(trees)
}
