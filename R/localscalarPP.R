#' localscalarPP
#'
#' A function that takes the output of a kappa, lambda, delta, VRates etc. RJ bayesTraits run and runs post-processing on it.
#' @param rjlog The RJ output of the run - typically suffixed with .VarRates.txt
#' @param tree The tree the analysis was run on
#' @param burnin The burnin (if required) for the mcmc (generally worked out from the other logfile)
#' @param thinning Thinning parameter for the MCMC output - again, worked out from the raw MCMC output logfile.
#' @export
#' @name localscalarPP

localscalarPP <- function(rjlog, rjtrees, tree, burnin = 0, thinning = 1) {
  # load the sample of trees.
  extree <- ladderize(tree)
  print("Loading log file.")
  rjout <- loadRJ(rjlog, burnin = burnin, thinning = thinning)
  print("Loading posterior trees.")
  posttrees <- read.nexus(rjtrees)
  posttrees <- posttrees[burnin:length(posttrees)]

  # Make a list to store descriptions of each scalar present in each iteration.
  ratesperit <- vector(mode = "list", length = nrow(rjout$rj_output))

  for (i in 1:nrow(rjout$rj_output)) {

    lastrates <- rjout$rj_output[i, !is.na(rjout$rj_output[i, ])]
    
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
  }

  # Get per-branch statistics.
  # Make a table called counts - this has one row per branch, the ancestor and descendant node, and then
  # space for the number of scalars in total, the type of scalars, and the number of generations one or more scalar
  # is applied, as well as space for those scalar types.
  counts <- matrix(ncol = 45, nrow = nrow(extree$edge))

  colnames(counts) <- c("branch", "ancNode", "descNode", "nTips", "start", "end", "mid", "orgBL", "newBL", "ratioBL", 
      "itersScaled", "itersRatescaled", "itersDelta", "itersKappa", "itersLambda", 
      "nScalar", "nRate", "nDelta", "nKappa", "nLambda",
      "rangeRate", "lqRate", "uqRate", "meanRate", "medianRate", "modeRate",
      "rangeDelta", "lqDelta", "uqDelta", "meanDelta", "medianDelta", "modeDelta",
      "rangeKappa", "lqKappa", "uqKappa", "meanKappa", "medianKappa", "modeKappa",
      "rangeLambda", "lqLambda", "uqLambda", "meanLambda", "medianLambda", "modeLambda", "species")

  counts[ , "branch"] <- c(1:nrow(extree$edge))
  counts[ , "ancNode"] <- extree$edge[ , 1]
  counts[ , "descNode"] <- extree$edge[ , 2]
  counts[ , "orgBL"] <- extree$edge.length
  print("Calculating mean branch lengths.")
  counts[ , "newBL"] <- meanBranches(reftree = extree, trees = posttrees)$meantree$edge.length
  counts[ , "ratioBL"] <- counts[ , "newBL"] / counts[ , "orgBL"]

  hts <- nodeHeights(extree)
  hts <- round(abs(hts - max(hts)), 4)
  counts[ , "start"] <- hts[ , 1]
  counts[ , "end"] <- hts[ , 2]
  counts[ , "mid"] <- mean(c(hts[ , 1]))
  counts <- as.data.frame(counts)

  for (i in 1:nrow(counts)) {
    descs <- getDescs(extree, node = counts[i, "descNode"])
    counts[i, "nTips"] <- sum(descs <= length(tree$tip.label))
    if (counts[i, "nTips"] == 0) {
      counts[i, "nTips"] <- 1
    }
    if (counts[i, "descNode"] <= length(extree$tip.label)) {
      counts[i, "species"] <- extree$tip.label[counts[i, "descNode"]]
    }
    counts[i, "mid"] <- mean(c(hts[i, 1], hts[i, 2]))
  }
  
  counts[ , c(11:44)] <- 0

  # Make tables to store the individual deltas and whatever for each iteration.
  # How to do this? I can't just have onc cell per branch per iteration - that's no good.
  # These ought to be lists I think, with one element per branch. Then add to each list each 
  # delta, or whatever. It will be slow to reasign, but I don't know how long it needs to be - potentially
  # very long... Especially big trees, long runs etc.

  rates <- vector(mode = "list", length = nrow(counts))
  deltas <- vector(mode = "list", length = nrow(counts))
  kappas <- vector(mode = "list", length = nrow(counts))
  lambdas <- vector(mode = "list", length = nrow(counts))

  for (i in 1:length(ratesperit)) {
    print(paste("Iteration", i, "of", length(ratesperit)))
    scalars <- ratesperit[[i]]
    # make a list of all the branches, as defined by descendent node, with a column of zeroes for any scalar, rate, delta, 
    # kappa and lambda. These become ones when a scalar is place.
     if (!all(is.na(ratesperit[[i]]))) {
      # Generate the table to store the current scalar of the current iteration.
      comptable <- matrix(ncol = 6, nrow = nrow(extree$edge))
      comptable[ , 1] <- extree$edge[ , 2]  
      comptable[ , c(2:6)] <- 0
      colnames(comptable) <- c("descNode", "totalscalar", "rate", "delta", "kappa", "lambda")
      
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
        mrca <- getMRCA(extree, rjout$taxa[rjout$taxa[ , 1] %in% taxa, 2])
        
        # Find out what type of scalar is currently under consideration.
        if (currenttrans == "Node") {
          # If it's a node, the descendents of that node, plus the branch leading to it, need to be considered.
          descs <- c(getDescs(extree, mrca), mrca)
          # Then find branches from descs and count.
          counts[counts[, "descNode"] %in% descs , "nScalar"] <- counts[counts[, "descNode"] %in% descs , "nScalar"] + 1
          counts[counts[, "descNode"] %in% descs , "nRate"] <- counts[counts[, "descNode"] %in% descs , "nRate"] + 1          
          rws <- comptable[ , 1] %in% descs
          comptable[rws & comptable[ , "totalscalar"] == 0, "totalscalar"] <- 1
          comptable[rws & comptable[ , "rate"] == 0, "rate"] <- 1
          # Finally record the rates.
          for (k in which(rws)) {
            rates[[k]] <- c(rates[[k]], currentscale)
          }
        }
        if (currenttrans == "Branch") {
          # Just use the branch that ends in mrca.
          descs <- mrca
          rws <- comptable[ , 1] %in% descs
          counts[counts[, "descNode"] == mrca , "nScalar"] <- counts[counts[, "descNode"] == mrca , "nScalar"] + 1
          counts[counts[, "descNode"] %in% descs , "nRate"] <- counts[counts[, "descNode"] %in% descs , "nRate"] + 1          
          comptable[rws & comptable[ , "totalscalar"] == 0, "totalscalar"] <- 1
          comptable[rws & comptable[ , "rate"] == 0, "rate"] <- 1
          # Finally record the deltas.
          for (k in which(rws)) {
            rates[[k]] <- c(rates[[k]], currentscale)
          }
        }
        if (currenttrans == "Delta") {
          # Define the branches (as descendent nodes) that were scaled.
          descs <- getDescs(extree, mrca)
          # Get rows of the table that need to have the sount increased.
          rws <- comptable[ , 1] %in% descs
          # Increase the number of deltas for each node in this scalar by 1, in the counts table, and 
          # the total number of scalars.
          counts[counts[, "descNode"] %in% descs , "nScalar"] <- counts[counts[, "descNode"] %in% descs , "nScalar"] + 1
          counts[counts[, "descNode"] %in% descs , "nDelta"] <- counts[counts[, "descNode"] %in% descs , "nDelta"] + 1
          # Then adjust the comptable to a) record yes or no for being scaled at all, and b) record 
          # yes or no for each of the scalars.
          comptable[rws & comptable[ , "totalscalar"] == 0, "totalscalar"] <- 1
          comptable[rws & comptable[ , "delta"] == 0, "delta"] <- 1

          # Finally record the deltas.
          for (k in which(rws)) {
            deltas[[k]] <- c(deltas[[k]], currentscale)
          }
        }
        if (currenttrans == "Kappa") {
          descs <- getDescs(extree, mrca)
          rws <- comptable[ , 1] %in% descs
          counts[counts[, "descNode"] %in% descs , "nScalar"] <- counts[counts[, "descNode"] %in% descs , "nScalar"] + 1
          counts[counts[, "descNode"] %in% descs , "nKappa"] <- counts[counts[, "descNode"] %in% descs , "nKappa"] + 1
          comptable[rws & comptable[ , "totalscalar"] == 0, "totalscalar"] <- 1
          comptable[rws & comptable[ , "kappa"] == 0, "kappa"] <- 1
          # Finally record the kappas.
          for (k in which(rws)) {
            kappas[[k]] <- c(kappas[[k]], currentscale)
          }
        }    
        if (currenttrans == "Lambda") {
          descs <- getDescs(extree, mrca)
          rws <- comptable[ , 1] %in% descs
          counts[counts[, "descNode"] %in% descs , "nScalar"] <- counts[counts[, "descNode"] %in% descs , "nScalar"] + 1
          counts[counts[ , "descNode"] %in% descs , "nLambda"] <- counts[counts[, "descNode"] %in% descs , "nLambda"] + 1
          comptable[rws & comptable[ , "totalscalar"] == 0, "totalscalar"] <- 1
          comptable[rws & comptable[ , "lambda"] == 0, "lambda"] <- 1
          # Finally record the lambdas.
          for (k in which(rws)) {
            lambdas[[k]] <- c(lambdas[[k]], currentscale)
          }          
        } 
      }      

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
  } 
  # Finally generate the descriptive stuff for the scalar values, and remove columns that are all zero.
  # First add in 1 where there is a NULL.

  for (i in 1:length(rates)) {
    if (is.null(rates[[i]])) {
      rates[[i]] <- 1
    }     
      counts[i, "rangeRate"] <- max(rates[[i]]) - min(rates[[i]])
      counts[i, "lqRate"] <- quantile(rates[[i]])[2]
      counts[i, "uqRate"] <- quantile(rates[[i]])[4]
      counts[i, "meanRate"] <- mean(rates[[i]])
      counts[i, "medianRate"] <- median(rates[[i]])
      if (length(rates[[i]]) > 1) {
        dens <- density(rates[[i]])
        counts[i, "modeRate"] <- dens$x[which(dens$y == max(dens$y))]
      } else {
        counts[i, "modeRate"] <- NA
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
        counts[i, "modeDelta"] <- dens$x[which(dens$y == max(dens$y))]
      } else {
        counts[i, "modeDelta"] <- NA
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
        counts[i, "modeKappa"] <- dens$x[which(dens$y == max(dens$y))]
      } else {
        counts[i, "modeKappa"] <- NA
      }
  }

  for (i in 1:length(lambdas)) {
    if (!is.null(lambdas[[i]])) {
      lambdas[[i]] <- 1
    }
      counts[i, "rangeLambda"] <- max(lambdas[[i]]) - min(lambdas[[i]])
      counts[i, "lqLambda"] <- quantile(lambdas[[i]])[2]
      counts[i, "uqLambda"] <- quantile(lambdas[[i]])[4]
      counts[i, "meanLambda"] <- mean(lambdas[[i]])
      counts[i, "medianLambda"] <- median(lambdas[[i]])
      if (length(lambdas[[i]]) > 1) {  
        dens <- density(lambdas[[i]])
        counts[i, "modeLambda"] <- dens$x[which(dens$y == max(dens$y))]
      } else {
        counts[i, "modeLambda"] <- NA
      }
  }

  # Finally remove zero columns.
  keeps <- NULL
  for (i in 1:ncol(counts)) {
    if (all(is.na(counts[ , i]))) {
      counts[ , i] <- 0
    }
  }

  for (i in 1:ncol(counts)) {
    if (any(counts[ , i] != 0)) {
      keeps <- c(keeps, i)
    }
  }
  counts <- counts[ , keeps]
  return(counts)
}
