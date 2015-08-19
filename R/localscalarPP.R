#' localscalarPP
#'
#' A function that takes the output of a kappa, lambda, delta, VRates etc. RJ bayesTraits run and runs post-processing on it.
#' @param rjlog The RJ output of the run - typically suffixed with .VarRates.txt
#' @param tree The tree the analysis was run on
#' @param burnin The burnin (if required) for the mcmc (generally worked out from the other logfile)
#' @param thinning Thinning parameter for the MCMC output - again, worked out from the raw MCMC output logfile.
#' @export
#' @name localscalarPP

localscalarPP <- function(rjlog, tree, burnin = 0, thinning = 1) {
  # load the sample of trees.
  extree <- tree
  rjout <- loadRJ(rjlog, burnin = burnin, thinning = thinning)

  ratesperit <- vector(mode = "list", length = nrow(rjout$rj_output))

  for (i in 1:nrow(rjout$rj_output)) {

    lastrates <- rjout$rj_output[i, !is.na(rjout$rj_output[i, ])]
    
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
  counts <- matrix(ncol = 13, nrow = nrow(extree$edge))

  counts[ , 1] <- c(1:nrow(extree$edge))
  counts[ , 2] <- extree$edge[ , 1]
  counts[ , 3] <- extree$edge[ , 2]
  counts[ , c(4:13)] <- 0

  colnames(counts) <- c("branch", "ancNode", "descNode", "timesScaled", "timesRatescaled", "timesDelta", "timesKappa", "timesLambda", 
    "nScalar", "nRate", "nDelta", "nKappa", "nLambda")

  for (i in 1:length(ratesperit)) {
    rates <- ratesperit[[i]]

    # make a list of all the branches, as defined by descendent node, with a column of zeroes for any scalar, rate, delta, 
    # kappa and lambda. These become ones when a scalar is place.
    comptable <- matrix(ncol = 6, nrow = nrow(extree$edge))
    comptable[ , 1] <- extree$edge[ , 2]  
    comptable[ , c(2:6)] <- 0
    colnames(comptable) <- c("descNode", "totalscalar", "rate", "delta", "kappa", "lambda")
    
    for (j in 1:nrow(rates)) {
      currentnode <- rates[j, "node"]
      currenttrans <- rates[j, "nodebranchdelta"]
      
      taxa <- rjout$subtrees[rjout$subtrees$node == currentnode, ]
      taxa <- taxa[ , !is.na(taxa)]
      taxa <- taxa[c(4:length(taxa))]
      
      mrca <- getMRCA(extree, rjout$taxa[rjout$taxa[ , 1] %in% taxa, 2])
      
      if (currenttrans == "Node") {
        descs <- c(getDescs(extree, mrca), mrca)
        # Then find branches from descs and count.
        counts[counts[, "descNode"] %in% descs , "nScalar"] <- counts[counts[, "descNode"] %in% descs , "nScalar"] + 1
        rws <- comptable[ , 1] %in% descs
        comptable[rws & comptable[ , "totalscalar"] == 0, "totalscalar"] <- 1
        comptable[rws & comptable[ , "rate"] == 0, "rate"] <- 1
      }
      if (currenttrans == "Branch") {
        # Just use the branch that ends in mrca.
        descs <- mrca
        rws <- comptable[ , 1] %in% descs
        counts[counts[, "descNode"] == mrca , "nScalar"] <- counts[counts[, "descNode"] == mrca , "nScalar"] + 1
        comptable[rws & comptable[ , "totalscalar"] == 0, "totalscalar"] <- 1
        comptable[rws & comptable[ , "rate"] == 0, "rate"] <- 1
      }
      if (currenttrans == "Delta") {
        descs <- getDescs(extree, mrca)
        rws <- comptable[ , 1] %in% descs
        counts[counts[, "descNode"] %in% descs , "nDelta"] <- counts[counts[, "descNode"] %in% descs , "nDelta"] + 1
        comptable[rws & comptable[ , "totalscalar"] == 0, "totalscalar"] <- 1
        comptable[rws & comptable[ , "delta"] == 0, "delta"] <- 1
      }    
      if (currenttrans == "Kappa") {
        descs <- getDescs(extree, mrca)
        rws <- comptable[ , 1] %in% descs
        counts[counts[, "descNode"] %in% descs , "nKappa"] <- counts[counts[, "descNode"] %in% descs , "nKappa"] + 1
        comptable[rws & comptable[ , "totalscalar"] == 0, "totalscalar"] <- 1
        comptable[rws & comptable[ , "kappa"] == 0, "kappa"] <- 1
      }    
      if (currenttrans == "Lambda") {
        descs <- getDescs(extree, mrca)
        rws <- comptable[ , 1] %in% descs
        counts[counts[, "descNode"] %in% descs , "nLambda"] <- counts[counts[, "descNode"] %in% descs , "nLambda"] + 1
        comptable[rws & comptable[ , "totalscalar"] == 0, "totalscalar"] <- 1
        comptable[rws & comptable[ , "lambda"] == 0, "lambda"] <- 1
      }
    
    }

    # Count up generations when scaling occurred.
    tmp <- counts[ , "descNode"] %in% comptable[comptable[ , "totalscalar"] == 1]
    counts[tmp, "timesScaled"] <- counts[tmp, "timesScaled"] + 1
    
    tmp <- counts[ , "descNode"] %in% comptable[comptable[ , "rate"] == 1]
    counts[tmp, "timesRatescaled"] <- counts[tmp, "timesRatescaled"] + 1
    
    tmp <- counts[ , "descNode"] %in% comptable[comptable[ , "delta"] == 1]
    counts[tmp, "timesDelta"] <- counts[tmp, "timesDelta"] + 1
    
    tmp <- counts[ , "descNode"] %in% comptable[comptable[ , "kappa"] == 1]
    counts[tmp, "timesKappa"] <- counts[tmp, "timesKappa"] + 1
        
    tmp <- counts[ , "descNode"] %in% comptable[comptable[ , "lambda"] == 1]
    counts[tmp, "timesLambda"] <- counts[tmp, "timesLambda"] + 1

  }
  return(counts)
}


### To add - 
# 1) 


# So this gives the number of times scaled, overall, per rate, and per transformation.
# The mean rate scalar needs to be added - this should be calculated per branch per generation (in the case of two scalars) and then over all per branch.
# A complication is that a rate scalar, then a delta transformation, is not really measureable or interpretable.
# Perhaps a measure of whether the deltas are >1 or <1 is useful as well - although again, in conjunction with rates I am not sure this is useful either...

# make a function to compare the branch lengths from the posterior to the original branch lengths as well - this could be done from the trees, or 
# from the raw output.

# Also there must be a way in which I can use the partition tree of the rj output, and then turn it into per-branch in a ape-style tree, this might be faster.








