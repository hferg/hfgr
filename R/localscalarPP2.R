#' localscalarPP2
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
#' @name localscalarPP2

localscalarPP2 <- function(rjlog, rjtrees, tree, burnin = 0, thinning = 1, returnscales = TRUE,
  returnorigins = TRUE) {

  extree <- ladderize(tree)
  print("Loading log file.")
  rjout <- loadRJ(rjlog, burnin = burnin, thinning = thinning)
    rj_output <- rjout$rj_output
    subtrees <- rjout$subtrees
    rjtaxa <- rjout$taxa
    niter <- nrow(rj_output)

  print("Loading posterior trees.")
  posttrees <- read.nexus(rjtrees)
  posttrees <- posttrees[burnin:length(posttrees)]
  print("Calculating mean branch lengths.")
  meanbl <- meanBranches(reftree = extree, trees = rjtrees, burnin = burnin, 
    thinning = thinning, pbar = TRUE)

  counts <- createCountsTable(extree, meanbl)

  # Make a list to store descriptions of each scalar present in each iteration.
  alltypes <- vector(mode = "list", length = nrow(rj_output))
  allmrcas <- vector(mode = "list", length = nrow(rj_output))

  rates <- matrix(rep(1, nrow(counts) * nrow(rj_output)), ncol = nrow(rj_output))
  rownames(rates) <- counts[ , "descNode"]

  # make lists for the origins of deltas etc.
  .tmp <- rep(1, nrow(rj_output))
  Node <- replicate(nrow(counts), as.numeric(paste(.tmp)), simplify = FALSE)
  Branch <- replicate(nrow(counts), as.numeric(paste(.tmp)), simplify = FALSE)
  Delta <- replicate(nrow(counts), as.numeric(paste(.tmp)), simplify = FALSE)
  Lambda <- replicate(nrow(counts), as.numeric(paste(.tmp)), simplify = FALSE)
  Kappa <- replicate(nrow(counts), as.numeric(paste(.tmp)), simplify = FALSE)
  names(Node) <- counts[ , "descNode"]
  names(Branch) <- counts[ , "descNode"]
  names(Delta) <- counts[ , "descNode"]
  names(Lambda) <- counts[ , "descNode"]
  names(Kappa) <- counts[ , "descNode"]

  print("Searching for scalars...")
  pb <- txtProgressBar(min = 0, max = nrow(rj_output), style = 3)
  for (i in 1:nrow(rj_output)) {
    lastrates <- rj_output[i, !is.na(rj_output[i, ])]
    
    # If the number of columns is seven, there are no scalars applied this generation.
    if (ncol(lastrates) == 7) {
      nodes[[i]] <- NA
      scales[[i]] <- NA
      types[[i]] <- NA
    } else {
      
      int <- lastrates[8:length(lastrates)]

      
      # TEST 
      nodes <- unlist(c(int[grep("NodeID*", names(int))]))
      scales <- unlist(c(int[grep("Scale*", names(int))]))
      types <- unlist(c(int[grep("NodeBranch*", names(int))]))
      taxa <- lapply(nodes, getTaxa)
      mrcas <- unlist(lapply(taxa, getMRCAhfg))
      alltypes[[i]] <- types
      allmrcas[[i]] <- mrcas
      # /TEST

      for (j in 1:length(mrcas)) {
        nm <- paste0(types[j], "[[\"", as.character(mrcas[j]), "\"]]", "[", i, "]")
        eval(parse(text = paste0(nm, "<-", scales[j])))
      }

      #scalars <- data.frame(node = nodes[[i]], scale = scales[[i]], type = types[[i]], mrca = unlist(mrcas[[i]]))
      
    }
    #origins <- fillOrigins(scalars = scalars, i = i, origins = origins)
    setTxtProgressBar(pb, i)    
  }
  origins <- list(nodes = Node, branch = Branch, delta = Delta,
    lambda = Lambda, kappa = Kappa, rates = rates)
  close(pb)

  # Make a list the same length for the taxa descendent from each node in nodes, and the mrca 
  # on the tree for each of those nodes. 

  alltypes <- unlist(alltypes)
  allmrcas <- unlist(allmrcas)

  # nOrgnScalar (i.e. origins of ANY scalar) can go now - since it doesn't matter what those scalars
  # are.

  bs <- table(unlist(allmrcas)[alltypes == "Branch"])
  ns <- table(unlist(allmrcas)[alltypes == "Node"])
  ds <- table(unlist(allmrcas)[alltypes == "Delta"])
  ks <- table(unlist(allmrcas)[alltypes == "Kappa"])
  ls <- table(unlist(allmrcas)[alltypes == "Lambda"])

  # count scalar origins.
  # This needs a little work to change the order of the bs, ns etc. to match the order in the counts table.

  counts$nOrgnBRate[counts$descNode %in% names(bs)] <- bs[match(counts$descNode[counts$descNode %in% names(bs)], names(bs))] 
  counts$nOrgnScalar[counts$descNode %in% names(bs)] <- counts$nOrgnBRate[counts$descNode %in% names(bs)] + 
    bs[match(counts$descNode[counts$descNode %in% names(bs)], names(bs))] 
  counts$nOrgnNRate[counts$descNode %in% names(ns)] <- ns[match(counts$descNode[counts$descNode %in% names(ns)], names(ns))]
  counts$nOrgnScalar[counts$descNode %in% names(ns)] <- counts$nOrgnBRate[counts$descNode %in% names(ns)] + 
    ns[match(counts$descNode[counts$descNode %in% names(ns)], names(ns))]
  counts$nOrgnDelta[counts$descNode %in% names(ds)] <- ds[match(counts$descNode[counts$descNode %in% names(ds)], names(ds))]
  counts$nOrgnScalar[counts$descNode %in% names(ds)] <- counts$nOrgnBRate[counts$descNode %in% names(ds)] + 
    ds[match(counts$descNode[counts$descNode %in% names(ds)], names(ds))]
  counts$nOrgnKappa[counts$descNode %in% names(ks)] <- ks[match(counts$descNode[counts$descNode %in% names(ks)], names(ks))]
  counts$nOrgnScalar[counts$descNode %in% names(ks)] <- counts$nOrgnBRate[counts$descNode %in% names(ks)] + 
    ks[match(counts$descNode[counts$descNode %in% names(ks)], names(ks))]
  counts$nOrgnLambda[counts$descNode %in% names(ls)] <- ls[match(counts$descNode[counts$descNode %in% names(ls)], names(ls))]
  counts$nOrgnScalar[counts$descNode %in% names(ls)] <- counts$nOrgnBRate[counts$descNode %in% names(ls)] + 
    ls[match(counts$descNode[counts$descNode %in% names(ls)], names(ls))]
  
  # Fill in transformation detail.

  counts[ , "meanDelta"] <- unlist(lapply(origins$delta, mean))
  counts[ , "medianDelta"] <- unlist(lapply(origins$delta, median))
  counts[ , "modeDelta"] <- unlist(lapply(origins$delta, modeStat))
  counts[ , "rangeDelta"] <- suppressWarnings(unlist(lapply(origins$delta, max)) - unlist(lapply(origins$delta, min)))
  counts[ , "meanKappa"] <- unlist(lapply(origins$kappa, mean))
  counts[ , "medianKappa"] <- unlist(lapply(origins$kappa, median))
  counts[ , "modeKappa"] <- unlist(lapply(origins$kappa, modeStat))
  counts[ , "rangeKappa"] <- suppressWarnings(unlist(lapply(origins$kappa, max)) - unlist(lapply(origins$kappa, min)))
  counts[ , "meanLambda"] <- unlist(lapply(origins$lambda, mean))
  counts[ , "medianLambda"] <- unlist(lapply(origins$lambda, median))
  counts[ , "modeLambda"] <- unlist(lapply(origins$lambda, modeStat))
  counts[ , "rangeLambda"] <- suppressWarnings(unlist(lapply(origins$lambda, max)) - unlist(lapply(origins$lambda, min)))
  counts[ , "meanRate"] <- apply(origins$rates, 1, mean)
  counts[ , "medianRate"] <- apply(origins$rates, 1, median)
  counts[ , "modeRate"] <- apply(origins$rates, 1, modeStat)
  counts[ , "rangeRate"] <- suppressWarnings(apply(origins$rates, 1, max) - unlist(apply(origins$rates, 1, min)))

  counts[ , "itersScaled"] <- 
  counts[ , "itersRatescaled"] <- apply(origins$rates, 1, function(x) sum(x != 1))
  counts[ , "itersDelta"] <- counts[ , "nOrgnDelta"]
  counts[ , "itersKappa"] <- counts[ , "nOrgnDelta"]
  counts[ , "itersLambda"] <- counts[ , "nOrgnDelta"]

  counts[ , "pScaled"] <- 
  counts[ , "pRate"] <- apply(origins$rates, 1, function(x) sum(x != 1)) / niter
  counts[ , "pDelta"] <- counts[ , "nOrgnDelta"] / niter
  counts[ , "pKappa"] <- counts[ , "nOrgnKappa"] / niter
  counts[ , "pLambda"] <- counts[ , "nOrgnLambda"] / niter

  counts[ , "nScalar"] <- 
  counts[ , "nRate"] <- counts[ , "nOrgnNRate"] + counts[ , "nOrgnBRate"]
  counts[ , "nDelta"] <- counts[ , "nOrgnDelta"]
  counts[ , "nKappa"] <- counts[ , "nOrgnDelta"]
  counts[ , "nLambda"] <- counts[ , "nOrgnDelta"]

  # Now just remove anything that is irrelevant (i.e. all values == 1) and then work out
  # what to return.

  counts <- counts[ , apply(counts, 2, function(x) all(x != 1))]
  counts <- counts[ , apply(counts, 2, function(x) all(x != 0))]

  meantree <- extree
  meantree$edge.length <- counts[c(2:nrow(counts)) , "meanBL"]

  res <- list(data = counts, niter = niter, meantree = meantree)

  if (sum(unlist((lapply(origins$delta, function(x) all(x != 1))))) == 0) {
    origins$delta <- NULL
  }

  if (sum(unlist((lapply(origins$kappa, function(x) all(x != 1))))) == 0) {
    origins$kappa <- NULL
  }

  if (sum(unlist((lapply(origins$lambda, function(x) all(x != 1))))) == 0) {
    origins$lambda <- NULL
  }

  if (sum(unlist((lapply(origins$nodes, function(x) all(x != 1))))) == 0) {
    origins$nodes <- NULL
  }  

  if (sum(unlist((lapply(origins$branches, function(x) all(x != 1))))) == 0) {
    origins$branches <- NULL
  }

  if (returnscales) {
    scalars <- list(rates = origins$rates)
    res <- c(res, list(scalars = scalars))
  }

  if (returnorigins) {
    res <- c(res, list(origins = origins))
  }

  return(res)
}
