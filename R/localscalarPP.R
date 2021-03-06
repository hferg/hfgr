#' localscalarPP
#'
#' A function that takes the output of a kappa, lambda, delta, VRates etc. RJ bayesTraits run and runs post-processing on it.
#' @param rjlog The RJ output of the run - typically suffixed with .VarRates.txt
#' @param tree The tree the analysis was run on
#' @param burnin The burnin (if required) for the mcmc (generally worked out from the other logfile)
#' @param thinning Thinning parameter for the MCMC output - again, worked out from the raw MCMC output logfile.
#' @param returnscales Return full distributions of all scalar values that hit each node?
#' @param returnorigins Return a list of the values of each scalar that originates at each node per iteration? 
#' @import phytools
#' @export
#' @name localscalarPP

localscalarPP <- function(rjlog, rjtrees, tree, burnin = 0, thinning = 1, returnscales = FALSE,
  returnorigins = FALSE) {
  # load the sample of trees.
  extree <- ladderize(tree)
  print("Loading log file.")
  rjout <- loadRJ(rjlog, burnin = burnin, thinning = thinning)
    rj_output <- rjout$rj_output
    subtrees <- subtrees
    rjtaxa <- rjout$taxa
  print("Loading posterior trees.")
  posttrees <- read.nexus(rjtrees)
  posttrees <- posttrees[burnin:length(posttrees)]

  # Make a list to store descriptions of each scalar present in each iteration.
  ratesperit <- vector(mode = "list", length = nrow(rj_output))

  print("Searching for scalars...")
  pb <- txtProgressBar(min = 0, max = nrow(rj_output), style = 3)
  for (i in 1:nrow(rj_output)) {

    lastrates <- rj_output[i, !is.na(rj_output[i, ])]
    
    # If the number of columns is seven, there are no scalars applied this generation.
    if (ncol(lastrates) == 7) {
      node <- NA
      scale <- NA
      created <- NA
      nodebranchdelta <- NA
    } else {
      
      int <- lastrates[8:length(lastrates)]

      node <- unlist(c(int[grep("NodeID*", names(int))]))
      scale <- unlist(c(int[grep("Scale*", names(int))]))
      created <- unlist(c(int[grep("CreatIt*", names(int))]))
      nodebranchdelta <- unlist(c(int[grep("NodeBranch*", names(int))]))
    
    }
    scalars <- data.frame(node = node, scale = scale, created = created, nodebranchdelta = nodebranchdelta)
    
    ratesperit[[i]] <- scalars
    setTxtProgressBar(pb, i)    
  }
  close(pb)


  counts <- matrix(ncol = 68, nrow = (nrow(extree$edge) + 1))

  colnames(counts) <- c("branch", "ancNode", "descNode", "nTips", "start", "end", "mid", "orgBL", "meanBL", "medianBL", "quart25", "quart75", 
      "itersScaled", "itersRatescaled", "itersDelta", "itersKappa", "itersLambda", 
      "pScaled", "pRate", "pDelta", "pKappa", "pLambda",
      "nScalar", "nRate", "nDelta", "nKappa", "nLambda",
      "nOrgnScalar", "nOrgnNRate", "nOrgnBRate", "nOrgnDelta", "nOrgnKappa", "nOrgnLambda",
      "mnSpI", "mnRpI", "mnDpI", "mnKpI", "mnLpI",
      "mnSpE", "mnRpE", "mnDpE", "mnKpE", "mnLpE",
      "rangeRate", "lqRate", "uqRate", "meanRate", "medianRate", "modeRate",
      "rangeDelta", "lqDelta", "uqDelta", "meanDelta", "medianDelta", "modeDelta",
      "rangeKappa", "lqKappa", "uqKappa", "meanKappa", "medianKappa", "modeKappa",
      "rangeLambda", "lqLambda", "uqLambda", "meanLambda", "medianLambda", "modeLambda", "species")

  counts[ , "branch"] <- c(0:nrow(extree$edge))
  counts[ , "ancNode"] <- c(0, extree$edge[ , 1])
  counts[ , "descNode"] <- c((length(tree$tip.label) + 1), extree$edge[ , 2])
  counts[ , "orgBL"] <- c(0, extree$edge.length)
  print("Calculating mean branch lengths.")
  meanbl <- meanBranches(reftree = extree, trees = rjtrees, burnin = burnin, thinning = thinning, pbar = TRUE)
  counts[ , "meanBL"] <- c(0, meanbl$meanbranches)
  counts[ , "medianBL"] <- c(0, meanbl$medianbranches)
  counts[ , "quart25"] <- c(0, meanbl$quart25)
  counts[ , "quart75"] <- c(0, meanbl$quart75)

  hts <- nodeHeights(extree)
  hts <- round(abs(hts - max(hts)), 4)
  counts[ , "start"] <- c(0, hts[ , 1])
  counts[ , "end"] <- c(0, hts[ , 2])
  counts <- as.data.frame(counts)

  # Deal with the root
  descs <- getDescs(extree, node = counts[1, "descNode"])
  counts[1, "nTips"] <- sum(descs <= length(tree$tip.label))
  counts[1, "mid"] <- 0
  counts[1, "species"] <- paste0(extree$tip.label[order(extree$tip.label)], collapse = ",")

  for (i in 2:nrow(counts)) {
    descs <- getDescs(extree, node = counts[i, "descNode"])
    counts[i, "nTips"] <- sum(descs <= length(tree$tip.label))
    if (counts[i, "nTips"] == 0) {
      counts[i, "nTips"] <- 1
    }
    if (counts[i, "descNode"] <= length(extree$tip.label)) {
      counts[i, "species"] <- extree$tip.label[counts[i, "descNode"]]
    } else {
      tips <- getDescs(extree, counts[i, "descNode"])
      tips <- tips[tips <= length(extree$tip.label)]
      tips <- extree$tip.label[tips]
      counts[i, "species"] <- paste0(sort(tips), collapse = ",")
    }
    counts[i, "mid"] <- mean(c(hts[(i - 1), 1], hts[(i - 1), 2]))
  }
  
  counts[ , c(13:67)] <- 0

  rates <- matrix(nrow = nrow(counts), ncol = length(ratesperit))
    rownames(rates) <- counts[ , "branch"]
  deltas <- vector(mode = "list", length = nrow(counts))
  kappas <- vector(mode = "list", length = nrow(counts))
  lambdas <- vector(mode = "list", length = nrow(counts))
  names(rates) <- c("root", extree$edge[ , 2])
  names(deltas) <- c("root", extree$edge[ , 2])
  names(kappas) <- c("root", extree$edge[ , 2])
  names(lambdas) <- c("root", extree$edge[ , 2])

  # make lists for the origins of deltas etc.
  .tmp <- rep(1, length(ratesperit))
  node_origins <- replicate(nrow(counts), as.numeric(paste(.tmp)), simplify = FALSE)
  branch_origins <- replicate(nrow(counts), as.numeric(paste(.tmp)), simplify = FALSE)
  delta_origins <- replicate(nrow(counts), as.numeric(paste(.tmp)), simplify = FALSE)
  lambda_origins <- replicate(nrow(counts), as.numeric(paste(.tmp)), simplify = FALSE)
  kappa_origins <- replicate(nrow(counts), as.numeric(paste(.tmp)), simplify = FALSE)
  names(node_origins) <- c("root", extree$edge[ , 2])
  names(branch_origins) <- c("root", extree$edge[ , 2])
  names(delta_origins) <- c("root", extree$edge[ , 2])
  names(lambda_origins) <- c("root", extree$edge[ , 2])
  names(kappa_origins) <- c("root", extree$edge[ , 2])

  print("Processing scalars...")
  pb <- txtProgressBar(min = 0, max = length(ratesperit), style = 3)
  for (i in 1:length(ratesperit)) {
    scalars <- ratesperit[[i]]
    # make a list of all the branches, as defined by descendent node, with a column of zeroes for any scalar, rate, delta, 
    # kappa and lambda. These become ones when a scalar is place.
     if (!all(is.na(ratesperit[[i]]))) {
      # Generate the table to store the current scalar of the current iteration.
      comptable <- matrix(ncol = 6, nrow = (nrow(extree$edge) + 1))
      comptable[1, 1] <- length(tree$tip.label) + 1
      comptable[c(2:nrow(comptable)) , 1] <- extree$edge[ , 2]  
      comptable[ , c(2:6)] <- 0
      colnames(comptable) <- c("descNode", "totalscalar", "rate", "delta", "kappa", "lambda")
      
      iteration_rates <- rep(1, nrow(counts))

      for (j in 1:nrow(scalars)) {
        # Get the node the scalar applies to, and the type of scalar.
        currentnode <- scalars[j, "node"]      
        currenttrans <- scalars[j, "nodebranchdelta"]
        currentscale <- as.numeric(as.character(scalars[j, "scale"]))

        # Find the taxa numbers that descend from the node that the scalar is applied too.
        taxa <- rjout$subtrees[rjout$subtrees$node == currentnode, ]
        taxa <- taxa[ , !is.na(taxa)]
        taxa <- taxa[c(4:length(taxa))]
        
        # Find the MRCA of those taxa - now the node number (mrca) is in terms of a phylogeny as stored in ape.
        # If the taxa is a single tip, make that MRCA (this can only apply to branch scalars.)
        if (length(taxa) == 1) {
          mrca <- which(tree$tip.label == rjout$taxa[rjout$taxa[ , 1] %in% taxa, 2])
        } else {
          mrca <- getMRCA(extree, rjout$taxa[rjout$taxa[ , 1] %in% taxa, 2])
        }

        # Find out what type of scalar is currently under consideration.
        if (currenttrans == "Node") {
          # If it's a node, the descendents of that node, plus the branch leading to it, need to be considered.
          descs <- c(getDescs(extree, mrca), mrca)
          # Then find branches from descs and count.
          counts[counts[, "descNode"] %in% descs , "nScalar"] <- counts[counts[, "descNode"] %in% descs , "nScalar"] + 1
          counts[counts[, "descNode"] %in% descs , "nRate"] <- counts[counts[, "descNode"] %in% descs , "nRate"] + 1 
          counts[counts[, "descNode"] == mrca, "nOrgnScalar"] <- counts[counts[, "descNode"] == mrca, "nOrgnScalar"] + 1
          counts[counts[, "descNode"] == mrca, "nOrgnNRate"] <- counts[counts[, "descNode"] == mrca, "nOrgnNRate"] + 1
          rws <- comptable[ , 1] %in% descs
          comptable[rws & comptable[ , "totalscalar"] == 0, "totalscalar"] <- 1
          comptable[rws & comptable[ , "rate"] == 0, "rate"] <- 1
          # Finally record the rates.
          if (mrca == length(extree$tip.label) + 1) {
            node_origins$root[i] <- currentscale
          } else {
            node_origins[[which(names(node_origins) == mrca)]][i] <- currentscale
          }

          for (k in which(rws)) {
              iteration_rates[k] <- iteration_rates[[k]] * currentscale
          }
        }
        if (currenttrans == "Branch") {
          # Just use the branch that ends in mrca.
          descs <- mrca
          rws <- comptable[ , 1] %in% descs
          counts[counts[, "descNode"] == mrca , "nScalar"] <- counts[counts[, "descNode"] == mrca , "nScalar"] + 1
          counts[counts[, "descNode"] %in% descs , "nRate"] <- counts[counts[, "descNode"] %in% descs , "nRate"] + 1
          counts[counts[, "descNode"] == mrca, "nOrgnScalar"] <- counts[counts[, "descNode"] == mrca, "nOrgnScalar"] + 1
          counts[counts[, "descNode"] == mrca, "nOrgnBRate"] <- counts[counts[, "descNode"] == mrca, "nOrgnBRate"] + 1          
          comptable[rws & comptable[ , "totalscalar"] == 0, "totalscalar"] <- 1
          comptable[rws & comptable[ , "rate"] == 0, "rate"] <- 1
          # Finally record the rates.
          if (mrca == length(extree$tip.label) + 1) {
            branch_origins$root[i] <- currentscale
          } else {
            branch_origins[[which(names(branch_origins) == mrca)]][i] <- currentscale
          }          
          for (k in which(rws)) {
              iteration_rates[k] <- iteration_rates[[k]] * currentscale
          }
        }
        if (currenttrans == "Delta") {
          # Define the branches (as descendent nodes) that were scaled.
          descs <- c(getDescs(extree, mrca), mrca)
          # Get rows of the table that need to have the sount increased.
          rws <- comptable[ , 1] %in% descs
          # Increase the number of deltas for each node in this scalar by 1, in the counts table, and 
          # the total number of scalars.
          counts[counts[, "descNode"] %in% descs , "nScalar"] <- counts[counts[, "descNode"] %in% descs , "nScalar"] + 1
          counts[counts[, "descNode"] %in% descs , "nDelta"] <- counts[counts[, "descNode"] %in% descs , "nDelta"] + 1
          counts[counts[, "descNode"] == mrca, "nOrgnScalar"] <- counts[counts[, "descNode"] == mrca, "nOrgnScalar"] + 1
          counts[counts[, "descNode"] == mrca, "nOrgnDelta"] <- counts[counts[, "descNode"] == mrca, "nOrgnDelta"] + 1
          # Then adjust the comptable to a) record yes or no for being scaled at all, and b) record 
          # yes or no for each of the scalars.
          comptable[rws & comptable[ , "totalscalar"] == 0, "totalscalar"] <- 1
          comptable[rws & comptable[ , "delta"] == 0, "delta"] <- 1

          # Finally record the deltas.
          if (mrca == length(extree$tip.label) + 1) {
            delta_origins$root[i] <- currentscale
          } else {
            delta_origins[[which(names(delta_origins) == mrca)]][i] <- currentscale
          }
          for (k in which(rws)) {
            deltas[[k]] <- c(deltas[[k]], currentscale)
          }
        }
        if (currenttrans == "Kappa") {
          descs <- c(getDescs(extree, mrca), mrca)
          rws <- comptable[ , 1] %in% descs
          counts[counts[, "descNode"] %in% descs , "nScalar"] <- counts[counts[, "descNode"] %in% descs , "nScalar"] + 1
          counts[counts[, "descNode"] %in% descs , "nKappa"] <- counts[counts[, "descNode"] %in% descs , "nKappa"] + 1
          counts[counts[, "descNode"] == mrca, "nOrgnScalar"] <- counts[counts[, "descNode"] == mrca, "nOrgnScalar"] + 1
          counts[counts[, "descNode"] == mrca, "nOrgnKappa"] <- counts[counts[, "descNode"] == mrca, "nOrgnKappa"] + 1
          comptable[rws & comptable[ , "totalscalar"] == 0, "totalscalar"] <- 1
          comptable[rws & comptable[ , "kappa"] == 0, "kappa"] <- 1
          # Finally record the kappas.
          if (mrca == length(extree$tip.label) + 1) {
            kappa_origins$root[i] <- currentscale
          } else {
            kappa_origins[[which(names(kappa_origins) == mrca)]][i] <- currentscale
          }          
          for (k in which(rws)) {
            kappas[[k]] <- c(kappas[[k]], currentscale)
          }
        }    
        if (currenttrans == "Lambda") {comptable
          descs <- c(getDescs(extree, mrca), mrca)
          rws <- comptable[ , 1] %in% descs
          counts[counts[, "descNode"] %in% descs , "nScalar"] <- counts[counts[, "descNode"] %in% descs , "nScalar"] + 1
          counts[counts[ , "descNode"] %in% descs , "nLambda"] <- counts[counts[, "descNode"] %in% descs , "nLambda"] + 1
          counts[counts[, "descNode"] == mrca, "nOrgnScalar"] <- counts[counts[, "descNode"] == mrca, "nOrgnScalar"] + 1
          counts[counts[, "descNode"] == mrca, "nOrgnLambda"] <- counts[counts[, "descNode"] == mrca, "nOrgnLambda"] + 1
          comptable[rws & comptable[ , "totalscalar"] == 0, "totalscalar"] <- 1
          comptable[rws & comptable[ , "lambda"] == 0, "lambda"] <- 1
          # Finally record the lambdas.
          if (mrca == length(extree$tip.label) + 1) {
            lambda_origins$root[i] <- currentscale
          } else {
            lambda_origins[[which(names(lambda_origins) == mrca)]][i] <- currentscale
          }          
          for (k in which(rws)) {
            lambdas[[k]] <- c(lambdas[[k]], currentscale)
          }
        } 
      }
  
      rates[ , i] <- iteration_rates

      # Now put the yes/no data from comptable into counts in the times scaled/rate, delta etc. This counts the number of 
      # ITERATIONS when scaling occurred  .
      tmp <- counts[ , "descNode"] %in% comptable[comptable[ , "totalscalar"] == 1, "descNode"] 
      counts[tmp, "itersScaled"] <- counts[tmp, "itersScaled"] + 1
      
      tmp <- counts[ , "descNode"] %in% comptable[comptable[ , "rate"] == 1, "descNode"]
      counts[tmp, "itersRatescaled"] <- counts[tmp, "itersRatescaled"] + 1
      
      tmp <- counts[ , "descNode"] %in% comptable[comptable[ , "delta"] == 1, "descNode"]
      counts[tmp, "itersDelta"] <- counts[tmp, "itersDelta"] + 1
      
      tmp <- counts[ , "descNode"] %in% comptable[comptable[ , "kappa"] == 1, "descNode"]
      counts[tmp, "itersKappa"] <- counts[tmp, "itersKappa"] + 1
          
      tmp <- counts[ , "descNode"] %in% comptable[comptable[ , "lambda"] == 1, "descNode"]
      counts[tmp, "itersLambda"] <- counts[tmp, "itersLambda"] + 1
    }
    setTxtProgressBar(pb, i)    
  }
  close(pb)
  
  # Finally generate the descriptive stuff for the scalar values, and remove columns that are all zero.
  # Before all this it is pretty straight forward to calculate the probs for each scalar.
  counts[ , "pScaled"] <- counts[ , "itersScaled"] / length(ratesperit)
  counts[ , "pRate"] <- counts[ , "itersRatescaled"] / length(ratesperit)
  counts[ , "pDelta"] <- counts[ , "itersDelta"] / length(ratesperit)
  counts[ , "pKappa"] <- counts[ , "itersKappa"] / length(ratesperit)
  counts[ , "pLambda"] <- counts[ , "itersLambda"] / length(ratesperit)

  counts[ , "mnSpI"] <- counts[ , "nScalar"] / length(ratesperit)
  counts[ , "mnRpI"] <- counts[ , "nRate"] / length(ratesperit)
  counts[ , "mnDpI"] <- counts[ , "nDelta"] / length(ratesperit)
  counts[ , "mnKpI"] <- counts[ , "nKappa"] / length(ratesperit)
  counts[ , "mnLpI"] <- counts[ , "nLambda"] / length(ratesperit)

  counts[ , "mnSpE"] <- counts[ , "nScalar"] / counts[ , "itersScaled"]
  counts[ , "mnRpE"] <- counts[ , "nRate"] / counts[ , "itersRatescaled"]
  counts[ , "mnDpE"] <- counts[ , "nDelta"] / counts[ , "itersDelta"]
  counts[ , "mnKpE"] <- counts[ , "nKappa"] / counts[ , "itersKappa"]
  counts[ , "mnLpE"] <- counts[ , "nLambda"] / counts[ , "itersLambda"]


  # This now works through a table, not a list...
  for (i in 1:nrow(rates)) {
      counts[i, "rangeRate"] <- max(rates[i, ]) - min(rates[i, ])
      counts[i, "lqRate"] <- quantile(rates[i, ])[2]
      counts[i, "uqRate"] <- quantile(rates[i, ])[4]
      counts[i, "meanRate"] <- mean(rates[i, ])
      counts[i, "medianRate"] <- median(rates[i, ])
      if (length(unique(rates[i, ])) == 1){
        counts[i, "modeRate"] <- unique(rates[i, ])
      } else {
        dens <- density(rates[i, ])
        counts[i, "modeRate"] <- dens$x[which(dens$y == max(dens$y))]
      }
  }

  for (i in 1:length(deltas)) {
    if (is.null(deltas[[i]])) {
      deltas[[i]] <- 1
    }      
      counts[i, "rangeDelta"] <- max(deltas[[i]]) - min(deltas[[i]])
      counts[i, "lqDelta"] <- quantile(deltas[[i]])[2]
      counts[i, "uqDelta"] <- quantile(deltas[[i]])[4]
      counts[i, "meanDelta"] <- mean(deltas[[i]])
      counts[i, "medianDelta"] <- median(deltas[[i]])
      if (length(deltas[[i]]) > 1) {
        dens <- density(deltas[[i]])
        if (length(dens$x[which(dens$y == max(dens$y))]) != 1) {
          dens_poss <- dens$x[which(dens$y == max(dens$y))]
          counts[i, "modeDelta"] <- sample(dens_poss, 1)
        } else {
          counts[i, "modeDelta"] <- dens$x[which(dens$y == max(dens$y))]
        }
      } else {
        counts[i, "modeDelta"] <- deltas[[i]]
      }
  }

  for (i in 1:length(kappas)) {
    if (is.null(kappas[[i]])) {
      kappas[[i]] <- 1
    }
      counts[i, "rangeKappa"] <- max(kappas[[i]]) - min(kappas[[i]])
      counts[i, "lqKappa"] <- quantile(kappas[[i]])[2]
      counts[i, "uqKappa"] <- quantile(kappas[[i]])[4]
      counts[i, "meanKappa"] <- mean(kappas[[i]])
      counts[i, "medianKappa"] <- median(kappas[[i]])
      if (length(kappas[[i]]) > 1) {
        dens <- density(kappas[[i]])
        if (length(dens$x[which(dens$y == max(dens$y))]) != 1) {
          dens_poss <- dens$x[which(dens$y == max(dens$y))]
          counts[i, "modeKappa"] <- sample(dens_poss, 1)
        } else {
          counts[i, "modeKappa"] <- dens$x[which(dens$y == max(dens$y))]
        }
      } else {
        counts[i, "modeKappa"] <- kappas[[i]]
      }
  }

  for (i in 1:length(lambdas)) {
    if (is.null(lambdas[[i]])) {
      lambdas[[i]] <- 1
    }
      counts[i, "rangeLambda"] <- max(lambdas[[i]]) - min(lambdas[[i]])
      counts[i, "lqLambda"] <- quantile(lambdas[[i]])[2]
      counts[i, "uqLambda"] <- quantile(lambdas[[i]])[4]
      counts[i, "meanLambda"] <- mean(lambdas[[i]])
      counts[i, "medianLambda"] <- median(lambdas[[i]])
      if (length(lambdas[[i]]) > 1) {  
        dens <- density(lambdas[[i]])
        if (length(dens$x[which(dens$y == max(dens$y))]) != 1) {
          dens_poss <- dens$x[which(dens$y == max(dens$y))]
          counts[i, "modeLambda"] <- sample(dens_poss, 1)
        } else {
          counts[i, "modeLambda"] <- dens$x[which(dens$y == max(dens$y))]
        }
      } else {
        counts[i, "modeLambda"] <- lambdas[[i]]
      }
  }

  # Finally remove zero columns.
  keeps <- NULL
  for (i in 1:ncol(counts)) {
    if (all(is.na(counts[ , i]))) {
      counts[ , i] <- 0
    }
    if (all(counts[ , i] == 1)) {
      counts[ , i] <- 0
    }    
  }

  for (i in 1:ncol(counts)) {
    if (any(counts[ , i] != 0)) {
      keeps <- c(keeps, i)
    }
  }
  counts <- counts[ , keeps]
  meantree <- extree
  meantree$edge.length <- counts[c(2:nrow(counts)) , "meanBL"]
  
  res <- list(data = counts, meantree = meantree, niter = length(ratesperit))

  if (returnscales) {
    scalars <- list(rate = rates, delta = deltas, kappa = kappas, lambda = lambdas)
    res <- c(res, list(scalars = scalars))
  }

  if (returnorigins) {
    origins <- list(nodes = node_origins, branches = branch_origins, delta = delta_origins, kappa = kappa_origins, lambda = lambda_origins)
    res <- c(res, list(origins = origins))
  }

  return(res)
}
