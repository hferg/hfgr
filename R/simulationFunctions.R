##############################################################################
#'
#' branchChangeSim
#' Simulate the expected change along a branch of length branchtime.
#' @export
#' @name branchChangeSim

branchChangeSim <- function(branchtime = 100, direction = "non", startval = 0, var = 0.2, increment = 1) { 
  # Create the output array, starting with the starting value
  simout <- array(dim = c(1,2), data = c(0, startval))
  varsimout <- NULL

  # Identify the time intervals within which we wish to simulate
  timeinc <- seq(from = 0, to = branchtime, by = increment) # divide the time up into increments
  if(max(timeinc) < branchtime) {
    timeinc <- c(timeinc, branchtime) # add on however long is left at the end. 
  } 

  # Loop through each of these time increments and simulate
  for (inc in 2:length(timeinc)) {
    multiplier <- timeinc[inc] - timeinc[inc - 1] # how much time happens along this increment?
    varsim  <- rnorm(n = 1, mean = 0, sqrt(var)) * multiplier # draw a random value from the distribution defined by var 
    # multiply by the amount of time that has occurred
    # Get the simulated value...
    if(direction == "neg") {
      simval <- simout[inc-1,2] - abs(varsim) 
    } else if(direction == "pos") {
      simval <- simout[inc-1,2] + abs(varsim) 
    } else if(direction == "non") {
      simval <- simout[inc-1,2] + varsim
    }

    # Add it to the table
    simout <- rbind(simout,c(timeinc[inc], simval))
  }
  return(simout)
}

##############################################################################
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

##############################################################################
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
#' @param lead Include the leading branch to the node in the simulations?
#' @name simOverTrees
#' @export

simOverTrees <- function(tree, reftrees, nsim, node, sigsq, direction = "non", increment = NULL, 
  lead = FALSE) {
  
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
        directional = direction, increment = increment, lead = lead)
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

##############################################################################
#' simuDatEb
#'
#' Simulate data over a phylogeny, or sub section of a phylogeny, in order to reflect variable rates.
#' Based on adding a constant to the tip data branch by branch and weighting by branch length and does
#' not need to stretch the phylogeny. This is good for exponential rate decay, but pretty shoddy for delta
#' at the moment - what it CAN do is do a good approximation of the delta relationship, but t isn't functioning
#' as half life at the moment.
#' @param tree The tree that data is simulated on.
#' @param node The node at which to make the changes.
#' @param a Initial rate increase before decay. If using the power decay this is, essentially, delta.
#' @param t In "pwr" (delta) mode this is ROUGHLY the delta parameter. In "exp" earyl burst mode, it is the number of half lives for the decay to take place over. Essentially a high number means the rate decays rapidly, and a low one means it decays slowly. This is used rather than a rate parameter, since decay is relative to branch lengths, but half lives are not.
#' @param sig The standard deviation of the random component of the BM simulation. NOT sigma squared.
#' @param decay Either "exp" (exponential) or "pwr" (power). The shape of the decay from the initial rate.
#' @param dat If you want to use some data you already have and change it this way, put it here.
#' @export
simuDatEB <- function(tree, node, sig, a, t, decay = "exp", dat = NULL) {
  
  if (isDefined(dat)) {
    dat <- dat
  } else {
    dat <- rTraitCont(tree, model = "BM", sigma = sig)
  }
  
  transdat <- dat  

  descs <- c(getDescs(tree, node), node)
  internal <- descs[descs > length(tree$tip.label)]
  tips <- descs[descs < length(tree$tip.label)]
  edges <- tree$edge[tree$edge[ ,2] %in% c(internal, tips), ]

  heights <- nodeHeights(tree)
  heights <- heights[tree$edge[ ,2] %in% descs, ]
  age <- max(heights)

  transnode <- nodeheight(tree, node)
  transtime <- age - transnode

  if (decay == "exp") {
    hl <- transtime / t
    r <- log(2) / hl
    fun <- function(x) a * (exp(1) ^ (-r * x))
  } else if (decay == "pwr") {
    if (t >= 1) {
      fun <- function(x) x ^ (1 / t)
    } else if (t < 1) {
      fun <- function(x) x ^ -(1 / t)
    }

  }
  
  for (i in 1:nrow(edges)) {
    tps <- getDescs(tree, edges[i, 2])
    tps <- tps[tps < length(tree$tip.label)]
    num <- sample(2, 1)
    
    start <- transtime - (age - nodeheight(tree, edges[i, 1]))
    end <- transtime - (age - nodeheight(tree, edges[i, 2]))
    
    m <- integrate(fun, start, end)$value
    
    if (num == 1) {
      transdat[tps] <- transdat[tps] + (m * sig)
    } else {
      transdat[tps] <- transdat[tps] - (m * sig)
    }
    
  }

  res <- list(original = dat, transformed = transdat)
  return(res)
}

##############################################################################
#' simuDatLambda
#'
#' Simulate data over a phylogeny, or part of a phylogeny, according to the model of Munkemuller et al 2012 by simulating data with a Brownian Motion process, and then adding noise to reflect loss of signal.
#' @param tree A tree of class 'phylo'
#' @param A node that describes the section of the tree to be transformed (if whole tree, the root node).
#' @param sig The standard deviation for the Brownian Motion process.
#' @param w The weight parameter. 1 = strong signal, 0 = weak signal.
#' @param dat A specified dataset, if you want to simulate the loss of signal on an existing dataset, ot something.

simuDatLambda <- function(tree, node, sig, w, dat = NULL) {
  
  if (isDefined(dat)) {
    dat <- dat
  } else {
    dat <- rTraitCont(tree, model = "BM", sigma = sig)
  }
  
  transdat <- dat  
  
  descs <- c(getDescs(tree, node), node)
  tips <- descs[descs < length(tree$tip.label)]
  tips <- tree$tip.label[tips]
  
  dd <- transdat[names(transdat) %in% tips]
  dd <- w * dd + (1 - w) * sample(dd)

  transdat[names(transdat) %in% tips] <- dd
  
  res <- list(original = dat, transformed = transdat)
  return(res)
}

##############################################################################
#' simuDatRate
#'
#' Simulate data over a phylogeny, or sub section of a phylogeny, in order to reflect variable rates.
#' Based on adding a constant to the tip data branch by branch and weighting by branch length and does
#' not need to stretch the phylogeny.
#' @name simuDatRate
#' @param tree The tree that data is simulated on.
#' @param transtree If you want to simulate data using a transformed tree to derive branch lengths, then this is where that tree goes. Defaults to NULL.
#' @param node The node at which to make the changes.
#' @param a The constant that is weighted by branch length and added to the data at the tips. Once the relationship between scalar and this is known, this can be scalar
#' @param sig The standard deviation of the random component of the BM simulation. NOT sigma squared.
#' @param dat If you want to use some data you already have and change it this way, put it here.
#' @param method Either add values recursively by tail branch length (default), or alternative "path" - which adds values proportional to total path length.
#' @param standardise Logical, whether or not to z-standardise the data. Defaults to FALSE.

simuDatRate <- function(tree, transtree = NULL, node, a, sig, dat = NULL, method = "recursive", standardise = FALSE) {
  
  if (isDefined(dat)) {
    dat <- dat
  } else {
    dat <- rTraitCont(tree, model = "BM", sigma = sig)
  }
  
  if (standardise == TRUE) {
    dat <- scale(dat, center = TRUE, scale = TRUE)
  }
  
  transdat <- dat

  if (is.null(transtree)) {
    transtree <- tree
  }

  if (node == "root") {
    node <- length(tree$tip.label) + 1
  }

  descs <- c(getDescs(transtree, node), node)
  internal <- descs[descs > length(transtree$tip.label)]
  tips <- descs[descs < length(transtree$tip.label)]
  edges <- transtree$edge[transtree$edge[ ,2] %in% c(internal, tips), ]
  bls <- transtree$edge.length[transtree$edge[ ,2] %in% c(internal, tips)]
  
  if (method == "recursive") {
    
    for (i in 1:nrow(edges)) {
      tps <- getDescs(transtree, edges[i, 2])
      tps <- tps[tps < length(transtree$tip.label)]
      num <- sample(2, 1)
      
      if (num == 1) {
        transdat[tps] <- transdat[tps] + (bls[i] * (a * sig))
      } else {
        transdat[tps] <- transdat[tps] - (bls[i] * (a * sig))
      }
    }
    
  } else if (method == "path") {
  
    print("Not implemented")
  
  }

  res <- list(original = dat, transformed = transdat)
  return(res)
}

##############################################################################
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

simTreeClades <- function(treesize, mintax, maxtax, clades, lambda = 0.2, mu = 0.02, 
  complete = FALSE, rescale = FALSE) {
    candidates <- NULL
    att <- 1
    if (clades * maxtax >= treesize) {
      stop("treesize is to small for the number and size of clades required - increase treesize.")
    }
    while (length(candidates) < clades) {
      print(paste("Attempt", att))
      tree <- sim.bd.taxa(treesize, numbsim = 1, lambda = lambda, mu = mu, complete = complete)[[1]]
      tree <- ladderize(tree)
      if (rescale == TRUE) {
        tree$edge.length <- tree$edge.length / max(nodeHeights(tree))
      } else if (is.numeric(rescale)) {
        tree$edge.length <- tree$edge.length / ((sum(tree$edge.length) / nrow(tree$edge)) / rescale)
      }
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

##############################################################################
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

# TODO: This is a prime candidate for refactoring!

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
