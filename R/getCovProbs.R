#' getCovProbs
#' A function to return the probabilities of independent or dependent at each node 
#' after fitting the covarion discrete model.
#' probgraphs will make a pdf of the overlap in probabilities at each node.
#' @export
#' @name getCovProbs

getCovProbs <- function(log, tree, ignore50 = FALSE, probgraphs = FALSE) {
  
  raw <- readLines(log)
  model <- gsub(" ", "", raw[2])
  
  mrca_location <- grep("MRCA", raw)
  header <- raw[c(mrca_location[1]:(grep("\\bTree Information\\b", raw) - 1))]
  mrca_location <- grep("MRCA", header)  
  
  output <- do.call(rbind, strsplit(raw[grep("\\bIteration\\b", raw):length(raw)], "\t"))
  colnames(output) <- output[1, ]
  output <- output[c(2:nrow(output)), ]
  output <- data.frame(output, stringsAsFactors = FALSE)
 
  for (i in 1:ncol(output)) {
    if (colnames(output)[i] != "Model.string" && colnames(output)[i] != "Dep...InDep") {
      output[ ,i] <- as.numeric(output[ ,i])
    }    
  }  

  posterior <- output

  cols <- colnames(posterior)
  nodes <- cols[grep("^X", cols)]
  xs <- strsplit(nodes, "\\.")
  nms <- vector(mode = "character", length = length(xs))
  
  for (i in 1:length(xs)) {
    nms[i] <- xs[[i]][1]
  }
  
  node_names <- unique(nms)
  res <- matrix(nrow = length(node_names), ncol = 3)
  colnames(res) <- c("tree_node", "I", "D")
  rownames(res) <- node_names

  if (probgraphs == TRUE) {
    pdfname <- paste0(log, ".nodeprobs.pdf")
    pdf(file = pdfname)
  }

  for (i in 1:length(node_names)) {
    if (i == length(mrca_location)) {
      tmp_t <- header[c(mrca_location[i]:length(header))]
    } else {
      tmp_t <- header[c(mrca_location[i]:(mrca_location[i + 1] - 1))]
    }

    mrca_name <- tmp_t[1]
    mrca_taxa <- tmp_t[c(2:length(tmp_t))]
    tmp_taxa <- vector(mode = "character", length = length(mrca_taxa))

    for (j in 1:length(mrca_taxa)) {
      tmp_taxa[j] <- strsplit(mrca_taxa, "\t")[[j]][2]
    }

    tmp <- posterior[ , grep(node_names[i], cols)]
    i_cols <- grep(paste0(node_names[i], "...I"), colnames(tmp))
    d_cols <- grep(paste0(node_names[i], "...D"), colnames(tmp))

    res[i, "tree_node"] <- findMRCA(tree, tmp_taxa)

    if (ignore50 == FALSE) {
      i_tmp <- rowSums(tmp[ , i_cols])
      d_tmp <- rowSums(tmp[ , d_cols])
      res[i, "I"] <- mean(i_tmp)
      res[i, "D"] <- mean(d_tmp)
    }

    if (ignore50 == TRUE) {
      i_tmp <- rowSums(tmp[ , i_cols])
      d_tmp <- rowSums(tmp[ , d_cols])
      i_tmp <- i_tmp[i_tmp != 0.5]
      d_tmp <- d_tmp[d_tmp != 0.5]
      res[i, "I"] <- mean(i_tmp)
      res[i, "D"] <- mean(d_tmp)
    }

    if(probgraphs == TRUE) {
      bwidth <- 3.5 * sd(rowSums(tmp[ , d_cols])) * length(rowSums(tmp[ , d_cols])) ^ -(1/3)
      hist <- ggplot(data.frame(parameter = rowSums(tmp[ , d_cols])), aes(x = parameter)) +
        geom_histogram(color = "darkgray", fill = "firebrick2", binwidth = bwidth) +
        geom_vline(xintercept = 0.5) +
        ggtitle(paste0("node ", res[i, "tree_node"], " (p(independent))"))
      print(hist)
    }
  }

  if (probgraphs == TRUE) {
    dev.off()
    print("Node probabilities pdf in working folder")
  }

  return(res)
}
