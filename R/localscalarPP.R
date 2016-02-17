#' localscalarPP
#'
#' A function that takes the output of a kappa, lambda, delta, VRates etc. RJ bayesTraits run and runs post-processing on it.
#' @param rjlog The RJ output of the run - typically suffixed with .VarRates.txt
#' @param tree The tree the analysis was run on
#' @param burnin The burnin (if required) for the mcmc (generally worked out from the other logfile)
#' @param thinning Thinning parameter for the MCMC output - again, worked out from the raw MCMC output logfile.
#' @param returnscales A vector of descendant nodes describing branches for which to return full distributions of scalars for.
#' @import phytools
#' @export
#' @name localscalarPP

localscalarPP <- function(rjlog, rjtrees, tree, burnin = 0, thinning = 1, returnscales = FALSE) {
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

  # ADJUST HERE IN ORDER TO ADD MEAN AND MEDIAN BLS INSTEAD OF NEW BLS.
  # Get per-branch statistics.
  # Make a table called counts - this has one row per branch, the ancestor and descendant node, and then
  # space for the number of scalars in total, the type of scalars, and the number of generations one or more scalar
  # is applied, as well as space for those scalar types.

  # TODO: I need to add in something to accommodate the root here. It will be the partition in Andrews thing with all
  #   taxa.

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
  meanbl <- meanBranches(reftree = extree, trees = posttrees, burnin = burnin, thinning = thinning)
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

  for (i in 2:nrow(counts)) {
    descs <- getDescs(extree, node = counts[i, "descNode"])
    counts[i, "nTips"] <- sum(descs <= length(tree$tip.label))
    if (counts[i, "nTips"] == 0) {
      counts[i, "nTips"] <- 1
    }
    if (counts[i, "descNode"] <= length(extree$tip.label)) {
      counts[i, "species"] <- extree$tip.label[counts[i, "descNode"]]
    }
    counts[i, "mid"] <- mean(c(hts[(i - 1), 1], hts[(i - 1), 2]))
  }
  
  counts[ , c(13:67)] <- 0

  # Make tables to store the individual deltas and whatever for each iteration.
  # How to do this? I can't just have onc cell per branch per iteration - that's no good.
  # These ought to be lists I think, with one element per branch. Then add to each list each 
  # delta, or whatever. It will be slow to reasign, but I don't know how long it needs to be - potentially
  # very long... Especially big trees, long runs etc.

  rates <- vector(mode = "list", length = nrow(counts))
  deltas <- vector(mode = "list", length = nrow(counts))
  kappas <- vector(mode = "list", length = nrow(counts))
  lambdas <- vector(mode = "list", length = nrow(counts))

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
  if (returnscales == TRUE) {
    res <- list(data = counts, rates = rates, deltas = deltas, kappas = kappas, lambdas = lambdas)
  } else {
    res <- counts
  }
  return(res)
}
