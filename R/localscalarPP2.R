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

localscalarPP2 <- function(rjlog, rjtrees, tree, burnin = 0, thinning = 1) {

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
  Node_effects <- replicate(nrow(counts), as.numeric(paste(.tmp)), simplify = FALSE)
  names(Node) <- counts[ , "descNode"]
  names(Branch) <- counts[ , "descNode"]
  names(Delta) <- counts[ , "descNode"]
  names(Lambda) <- counts[ , "descNode"]
  names(Kappa) <- counts[ , "descNode"]
  names(Node_effects) <- counts[ , "descNode"]

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
 
      nodes <- unlist(c(int[grep("NodeID*", names(int))]))
      scales <- unlist(c(int[grep("Scale*", names(int))]))
      types <- unlist(c(int[grep("NodeBranch*", names(int))]))
      taxa <- lapply(nodes, getTaxa, subtrees = subtrees)
      mrcas <- unlist(lapply(taxa, getMRCAhfg, tree = extree))
      alltypes[[i]] <- types
      allmrcas[[i]] <- mrcas

      for (j in 1:length(mrcas)) {
        nm <- paste0(types[j], "[[\"", as.character(mrcas[j]), "\"]]", "[", i, "]")
        eval(parse(text = paste0(nm, "<-", scales[j])))
      }

      #scalars <- data.frame(node = nodes[[i]], scale = scales[[i]], type = types[[i]], mrca = unlist(mrcas[[i]]))
      
    }
    #origins <- fillOrigins(scalars = scalars, i = i, origins = origins)
    setTxtProgressBar(pb, i)    
  }

  close(pb)

  # If there are node scalars, I need to find all the branches that descend from each node
  # in the Node_effects table, and multiply the existing scalar by the node one for each
  # generation. Then multiply by branch scalars as well.

  for (i in 1:length(Node)) {
    .tmp <- multiplyNodes(Node[[i]], names(Node)[i], extree, Node_effects)
    Node_effects[names(.tmp)] <- .tmp
  }

  Node_effects <- lapply(1:length(Node_effects), function(x) Node_effects[[x]] * Branch[[x]])
  names(Node_effects) <- counts[ , "descNode"]
  
  origins <- list(nodes = Node, branches = Branch, delta = Delta,
    lambda = Lambda, kappa = Kappa, rates = Node_effects)

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

  bstaxa <- counts$descNode %in% names(bs)
  nstaxa <- counts$descNode %in% names(ns)
  dstaxa <- counts$descNode %in% names(ds)
  kstaxa <- counts$descNode %in% names(ks)
  lstaxa <- counts$descNode %in% names(ls)

  counts$nOrgnBRate[bstaxa] <- bs[match(counts$descNode[bstaxa], names(bs))] 
  counts$nOrgnScalar[bstaxa] <- counts$nOrgnBRate[bstaxa] + bs[match(counts$descNode[bstaxa], names(bs))] 
  
  counts$nOrgnNRate[nstaxa] <- ns[match(counts$descNode[nstaxa], names(ns))]
  counts$nOrgnScalar[nstaxa] <- counts$nOrgnBRate[nstaxa] + ns[match(counts$descNode[nstaxa], names(ns))]
  
  counts$nOrgnDelta[dstaxa] <- ds[match(counts$descNode[dstaxa], names(ds))]
  counts$nOrgnScalar[dstaxa] <- counts$nOrgnBRate[dstaxa] + ds[match(counts$descNode[dstaxa], names(ds))]
  
  counts$nOrgnKappa[kstaxa] <- ks[match(counts$descNode[kstaxa], names(ks))]
  counts$nOrgnScalar[kstaxa] <- counts$nOrgnBRate[kstaxa] + ks[match(counts$descNode[kstaxa], names(ks))]
  
  counts$nOrgnLambda[lstaxa] <- ls[match(counts$descNode[lstaxa], names(ls))]
  counts$nOrgnScalar[lstaxa] <- counts$nOrgnBRate[lstaxa] + ls[match(counts$descNode[lstaxa], names(ls))]
  
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
  counts[ , "meanRate"] <- sapply(origins$rates, mean)
  counts[ , "medianRate"] <- sapply(origins$rates, median)
  counts[ , "modeRate"] <- sapply(origins$rates, modeStat)
  counts[ , "rangeRate"] <- suppressWarnings(sapply(origins$rates, max) - unlist(sapply(origins$rates, min)))

  counts[ , "itersScaled"] <- 
  counts[ , "itersRatescaled"] <- sapply(origins$rates, function(x) sum(x != 1))
  counts[ , "itersDelta"] <- counts[ , "nOrgnDelta"]
  counts[ , "itersKappa"] <- counts[ , "nOrgnDelta"]
  counts[ , "itersLambda"] <- counts[ , "nOrgnDelta"]

  counts[ , "pScaled"] <- 
  counts[ , "pRate"] <- sapply(origins$rates, function(x) sum(x != 1)) / niter
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

  if (all(sapply(origins$delta, function(x) all(x == 1)))) {
    origins$delta <- NULL
  }

  if (all(sapply(origins$kappa, function(x) all(x == 1)))) {
    origins$kappa <- NULL
  }

  if (all(sapply(origins$lambda, function(x) all(x == 1)))) {
    origins$lambda <- NULL
  }

  if (all(sapply(origins$nodes, function(x) all(x == 1)))) {
    origins$nodes <- NULL
  }  

  if (all(sapply(origins$branches, function(x) all(x == 1)))) {
    origins$branches <- NULL
  }

  if (all(sapply(origins$rates, function(x) all(x == 1)))) {
    origins$rates <- NULL
  }

  res <- list(data = counts, niter = niter, meantree = meantree)

  if (!is.null(origins$rates)) {
    scalars <- list(rates = origins$rates)
    origins$rates <- NULL
    res <- c(res, list(scalars = scalars))
  }

  res <- c(res, list(origins = origins))

  return(res)
}
