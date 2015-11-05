#' readDiscMaps
#'
#' A function to read in and interpret the discrete/independent maps that result from model 11.

readDiscMaps <- function(logfile, tree, burnin = 0, thinning = 1, returnpies = FALSE) {

  tree <- ladderize(tree)
  raw <- readLines(logfile)
  # Then find the node IDs, which are numbered according to map
  # position and labelled with "HNode

  hnodes <- raw[grep("HNode", raw)]
  # now I need to split that into a list, or something...

  # Find out how many samples there are in the posterior, to make a result table.
  nsamples <- nrow(do.call(rbind, strsplit(raw[grep("\\bIteration\\b", raw):length(raw)], "\t")))

  # These are all different lengths - so what I need to do is move through them, and find 
  # out which node is which using findMRCA, and just make a table that way.

  preres <- matrix(ncol = length(hnodes), nrow = nsamples)
  cols <- vector(mode = "character", length = ncol(preres))

  for (i in 1:length(hnodes)) {
    current <- strsplit(hnodes[[i]], "\t")[[1]]
    tips <- current[4:length(current)]
    mrca <- getMRCA(tree, tips)
    if (length(tips) > 1) {
      cols[i] <- mrca
    } else {
      cols[i] <- tips
    }
  }

  colnames(preres) <- cols

  # And now I need to unlist the maps according to , and add them into the table of results.
  # First find where the iterations start/
  start <- grep("\\bIteration\\b", raw)
  map_pos <- which(strsplit(raw[start], "\t")[[1]] == "Map")

  for (i in 1:nrow(preres)) {
    preres[i, ] <- strsplit(strsplit(raw[start + i], "\t")[[1]][map_pos], ",")[[1]]
  }

  # So now each column is all the 1 and 0 associated with that node.
  # Generate some summaries.

  res <- matrix(nrow = ncol(preres), ncol = 2)
  rownames(res) <- colnames(preres)
  colnames(res) <- c(0, 1)
  
  for (i in 1:ncol(preres)) {
    res[i, 1] <- table(preres[ , i])[1]
    res[i, 2] <- table(preres[ , i])[2]
  }

  if (returnpies == TRUE) {
    nodes <- res[!is.na(as.numeric(as.character(rownames(res)))), ]
    nodes <- nodes / sum(nodes[1, ])
    tmp <- c(1, 0)
    nodes <- rbind(tmp, nodes)    
    res <- list(allcounts = res, pies = nodes)
  }

  return(res)
}
