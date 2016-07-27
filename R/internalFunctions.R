#' 
#' #' A bunch of internal functions for hfgr, mostly for the localscalar PP.
#' 

createCountsTable <- function(extree, meanbl) {
  counts <- matrix(ncol = 52, nrow = (nrow(extree$edge) + 1))
  colnames(counts) <- c("branch", "ancNode", "descNode", "nTips", "start", "end", "mid", "orgBL", "meanBL", "medianBL", "quart25", "quart75", 
      "itersScaled", "itersRatescaled", "itersDelta", "itersKappa", "itersLambda", 
      "pScaled", "pRate", "pDelta", "pKappa", "pLambda",
      "nScalar", "nRate", "nDelta", "nKappa", "nLambda",
      "nOrgnScalar", "nOrgnNRate", "nOrgnBRate", "nOrgnDelta", "nOrgnKappa", "nOrgnLambda",
      "rangeRate", "lqRate", "uqRate", "meanRate", "medianRate", "modeRate",
      "rangeDelta", "meanDelta", "medianDelta", "modeDelta",
      "rangeKappa", "meanKappa", "medianKappa", "modeKappa",
      "rangeLambda", "meanLambda", "medianLambda", "modeLambda", "species")

  counts[ , "branch"] <- c(0:nrow(extree$edge))
  counts[ , "ancNode"] <- c(0, extree$edge[ , 1])
  counts[ , "descNode"] <- c((length(tree$tip.label) + 1), extree$edge[ , 2])
  counts[ , "orgBL"] <- c(0, extree$edge.length)
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
  counts[1, "nTips"] <- sum(descs <= length(extree$tip.label))
  counts[1, "mid"] <- 0
  counts[1, "species"] <- paste0(extree$tip.label[order(extree$tip.label)], collapse = ",")

  for (i in 2:nrow(counts)) {
    descs <- getDescs(extree, node = counts[i, "descNode"])
    counts[i, "nTips"] <- sum(descs <= length(extree$tip.label))
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
  
  counts[ , c(13:51)] <- 0
  return(counts)
}

getTaxa <- function(x, subtrees) {
  taxa <- subtrees[subtrees$node == x, ]
  taxa <- taxa[ , !is.na(taxa)]
  taxa <- taxa[c(4:length(taxa))]
  return(as.numeric(unlist(taxa)))
}

getMRCAhfg <- function(x, tree, rjtaxa) {
  if (length(x) == 1) {
    mrca <- which(tree$tip.label == rjtaxa[rjtaxa[ , 1] %in% x, 2])
  } else {
    mrca <- getMRCA(tree, rjtaxa[rjtaxa[ , 1] %in% x, 2])
  }
  return(mrca)
}

multiplyNodes <- function(scales, name, tree, Node_effects) {
  # get descendents
  descs <- getDescs(tree, name)
  .tmp <- lapply(Node_effects[as.character(descs)], function(x) x * scales)
  return(.tmp)
}
