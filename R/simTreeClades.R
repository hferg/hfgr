#' simTreeClades
#' 
#' This function simulates a tree of a given size, and then identifies a given number of clades of a certain size.
#' Trees are simulated using the function sim.bd.taxa from TreeSim. the function can get into an infinite loop
#' if you aren't careful about the specifications, for example if the tree has 100 tips, a single clade of 95
#' taxa might take a long time to find.
#' @param treesize The number of tips the final tree should have.
#' @param mintax The minimum size of each of the clades to be identified.
#' @param maxtax The maximum size of each of the clades to be identified.
#' @param caldes The number of clades to be identified.
#' @param lambda The birthrate for the tree simulation process, defaults to 0.2
#' @param mu The deathreate for the tree simulation process, defaults to 0.02
#' @param complete If FALSE, returns just the extant taxa at the end of the simulation, if TRUE extinct taxa are retained.
#' @return A list with two elements, the first is the tree, and the second is a list of the taxa names for each clade of the correct size identified.
#' @export

simTreeClades <- function(treesize, mintax, maxtax, clades, lambda = 0.2, mu = 0.02, complete = FALSE) {
    candidates <- NULL
    att <- 1
    if (clades * maxtax >= treesize) {
      stop("treesize is to small for the number and size of clades required - increase treesize.")
    }
    while (length(candidates) < clades) {
      print(paste("Attempt", att))
      tree <- sim.bd.taxa(treesize, numbsim = 1, lambda = lambda, mu = mu, complete = complete)[[1]]
      tree <- ladderize(tree)
      tree$edge.length <- tree$edge.length / max(nodeHeights(tree))
      nodes <- matrix(nrow = nrow(tree$edge), ncol = 2)
      colnames(nodes) <- c("Node", "nTips")
      nodes[ , 1] <- tree$edge[ , 2]

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
      att <- att + 1
    }

    if (clades == 1) {
      change_nodes <- candidates[sample(length(candidates), 1)]
    } else {
      change_nodes <- sample(candidates, clades)
    }

    change_taxa <-  vector(mode = "list", length = length(change_nodes))
    for (k in 1:length(change_taxa)) {
      tips <- getDescs(tree, change_nodes[[k]])
      tips <- tips[tips <= length(tree$tip.label)]
      tips <- tree$tip.label[tips]
      change_taxa[[k]] <- tips
    }
    trees <- list(tree = tree, changed_taxa = change_taxa)
  return(trees)
}
