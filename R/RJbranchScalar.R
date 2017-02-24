#' RJbranchScalars
#' A function that takes the output of a kappa, lambda, delta model and produces a summary of
#' which branches were scaled, how much, and what value that scalar has, and a tree with the mean branch
#' lengths of the trees sampled during the analysis.
#' @export
#' @param rjlog The logfile from the RJ portion of the BT analysis.
#' @param rjtrees The filename of the posterior sample of trees from the analysis.
#' @param tree The tree on which the analysis was performed.

RJbranchScalar <- function(rjlog, rjtrees, tree) {

  out <- loadRJ(rjlog)
  posttrees <- read.nexus(rjtrees)
  taxa <- out$taxa
  subtrees <- out$subtrees
  rj <- out$rj_output
  nodes_scalars <- list(nodes = rj[ ,grep("Node.ID", colnames(rj))], scalars = rj[ ,grep("Scale_", colnames(rj))])

  resArray <- array(list(NULL), c(nrow(tree$edge), (nrow(nodes_scalars$nodes))))

  for (i in 1:nrow(nodes_scalars$nodes)) {
    nodes <- nodes_scalars$nodes[i , !is.na(nodes_scalars$nodes[i, ])]
    scalars <- nodes_scalars$scalars[i, !is.na(nodes_scalars$scalars[i, ])]

    if (length(nodes) == 0) {

    } else if (length(nodes) > 1) {
      for (j in 1:ncol(nodes)) {
        taxcols <- which(!is.na(subtrees[subtrees$node == nodes[1, j], ]))
        tax <- subtrees[subtrees$node == nodes[1, j], taxcols[c(4:length(taxcols))]]
        tips <- taxa[as.numeric(tax[1, ]), 2]
        node <- getMRCA(tree, tips)

        descs <- getDescs(tree, node)
        edges <- which(tree$edge[ ,2] %in% descs)
        edges <- c(edges, which(tree$edge[ ,2] == node))

        for (k in edges) {
          if (is.null(resArray[[k, i]])) {
            resArray[[k, i]] <- as.numeric(scalars[1, j])
          } else {
            resArray[[k, i]] <- c(resArray[[k, i]] <- as.numeric(scalars[1, j]))
          }
        }
      }
    } else {
      taxcols <- which(!is.na(subtrees[subtrees$node == nodes, ]))
      tax <- subtrees[subtrees$node == nodes, taxcols[c(4:length(taxcols))]]
      tips <- taxa[as.numeric(tax[1, ]), 2]
      node <- getMRCA(tree, tips)

      descs <- getDescs(tree, node)
      edges <- which(tree$edge[ ,2] %in% descs)
      edges <- c(edges, which(tree$edge[ ,2] == node))

      for (k in edges) {
        if (is.null(resArray[[k, i]])) {
          resArray[[k, i]] <- as.numeric(scalars)
        } else {
          resArray[[k, i]] <- c(resArray[[k, i]] <- as.numeric(scalars))
        }
      }
    }
  }

  result <- matrix(ncol = 11, nrow = nrow(tree$edge))
  colnames(result) <- c("anc", "dec", "freq", "prop", "mn", "mn_wp", "md", "md_wp", "mod", "mod_wp", "taxa")
  result[ ,"anc"] <- tree$edge[ ,1]
  result[ ,"dec"] <- tree$edge[ ,2]
  result[ ,"freq"] <- 0
  result[ ,"prop"] <- 0
  result[ ,"mn"] <- 0
  result[ ,"mn_wp"] <- 0
  result[ ,"md"] <- 0
  result[ ,"md_wp"] <- 0
  result[ ,"mod"] <- 0
  result[ ,"mod_wp"] <- 0

  for (i in 1:nrow(result)) {
    result[i, "freq"] <- length(unlist(resArray[i, ]))
    result[i, "prop"] <- length(unlist(resArray[i, ])) / ncol(resArray)
    result[i, "mn"] <- mean(unlist(resArray[i, ]))
    result[i, "mn_wp"] <- mean(c(unlist(resArray[i, ]), rep(1, (ncol(resArray) - length(unlist(resArray[i, ]))))))
    
    if (length(unlist(resArray[i, ])) <= 2) {
      result[i, "md"] <- NA
      result[i, "md_wp"] <- NA
      result[i, "mod"] <- NA
      result[i, "mod_wp"] <- NA
    } else {
      result[i ,"mod"] <- modeStat(unlist(resArray[i, ]))
      result[i ,"mod_wp"] <- modeStat(c(unlist(resArray[i, ]), rep(1, (ncol(resArray) - length(unlist(resArray[i, ]))))))
      result[i ,"md"] <- median(unlist(resArray[i, ]))
      result[i ,"md_wp"] <- median(c(unlist(resArray[i, ]), rep(1, (ncol(resArray) - length(unlist(resArray[i, ]))))))      
    }
    
  }

  result <- data.frame(result)
  
  tree$edge.length <- meanBranches(posttrees)
  
  res <- list(tree = tree, branches = result, full_out = resArray)
  return(res)
}

