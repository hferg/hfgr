#' summariseModelString
#'
#' Standardises the model string output from an RJ analysis, and returns a list of two dataframes.
#' One of these gives the frequency of each unique model and the parameter composition of it, and 
#' the other details the proportion of times each pair of parameters are locked to the same value.
#' This currently doesn't care whether the parameters are fixed to be the same, or whether they are
#' both zero - just that they are both the same.
#' @param logfile The name of the logfile resulting from an RJ BayesTraits model.
#' @export
#' @name summariseModelString

summariseModelString <- function(logfile) {
  # get model strings
  output <- btmcmc(logfile)
  modstr <- output[ ,"Model.string"]
  
  # convert to a matrix.
  modstr <- gsub("'", "", modstr)
  modstr <- strsplit(modstr, " ")
  modstr <- matrix(unlist(modstr), ncol = length(modstr[[1]]), byrow = TRUE)
  
  pars <- colnames(output)[c((which(colnames(output) =="Model.string") + 1) 
    : (which(colnames(output) =="Model.string") + ncol(modstr)))]
    
  zerotable <- matrix(nrow = ncol(modstr), ncol = 2)
  rownames(zerotable) <- pars
  colnames(zerotable) <- c("frequency", "prob")
  
  for (i in 1:nrow(zerotable)) {
    zerotable[i, 1] <- sum(modstr[ ,i] == "Z")
    zerotable[i, 2] <- sum(modstr[ ,i] == "Z") / nrow(modstr)
  }
  
  # work out IDs.
  mods <- list()
  full.mods <- list()
  
  for (i in 1:nrow(modstr)) {
    ids <- list()
    
    for (j in 1:ncol(modstr)) {
      ids[[j]] <- sum(2 ^ (which(modstr[i, ] == modstr[i, j])))
      
      if (modstr[i, j] == "Z") {
        ids[[j]] <- paste0(ids[[j]], "Z")
      }
      
    }
    
    mods[[i]] <- unlist(unique(ids))
    full.mods[[i]] <- unlist(ids)
  }
  
  mods.cf <- mods
  
  for (i in 1:length(mods.cf)) {
    mods.cf[[i]] <- paste(mods.cf[[i]], collapse = "-")
  }

  full.mods <- lapply(unique(full.mods), setNames, nm = pars)
  names(full.mods) <- unlist(unique(mods.cf))
  res <- matrix(nrow = length(unique(mods)), ncol = max(sapply(mods, length)) + 4)
  colnames(res) <- c("id", "freq", "prop", "Z", paste0("g", c(1:max(sapply(mods, length)))))


  res[ ,"id"] <- names(full.mods)
  res[ ,"freq"] <- table(unlist(mods.cf))
  res[ ,"prop"] <- table(unlist(mods.cf)) / length(mods)

  for (i in 1:nrow(res)) {
    mod.conf <-  full.mods[names(full.mods) == res[i, "id"]][[1]]
    # col correct fits models to the correct column and is adjusted if the Z column is used.
    col.correct <- 4

    for (j in 1:length(strsplit(res[i, "id"], "-")[[1]])) {
      piece <- strsplit(res[i, "id"], "-")[[1]][j]
      group <- which(mod.conf == piece)

      if (grepl("Z", piece)) {
        res[i, "Z"] <- paste(names(group), collapse = "+")
        col.correct <- 3
      } else {
        res[i, (col.correct + j)] <- paste(names(group), collapse = "+")
      }
    }
    
  }

  # sort res based on frequency
  res <- data.frame(res)
  res$freq <- as.numeric(as.character(res$freq))
  res <- res[base::order(res$freq, decreasing = TRUE), ]
  rownames(res) <- c(1:nrow(res))
  
  # The last column may be all NA (since the colnum calculation doesn't account for Z)
  # so drop columns until it is not.
  while (all(is.na(res[ , length(res)]))) {
    res <- res[ ,c(1:length(res)-1)]
  }
  

  pairwise <- matrix(ncol = length(pars), nrow = length(pars), "-")
  rownames(pairwise) <- colnames(pairwise) <- pars
  pairwise.f <- pairwise
  pairwise.p <- pairwise
  pairwise.f.nz <- pairwise
  pairwise.p.nz <- pairwise
  pairwise.f.jz <- pairwise
  pairwise.p.jz <- pairwise
  
  pairs <- combn(pars, 2)
  
  for (i in 1:ncol(pairs)) {
    tot <- 0
    totnz <- 0
    totz <- 0

    for (j in 1:length(full.mods)) {
      pars.tmp <- which(names(full.mods[[j]]) %in% pairs[ ,i])

      if (full.mods[[j]][pars.tmp[1]] == full.mods[[j]][pars.tmp[2]]) {
      
        tot <- tot + res$freq[which(res$id == names(full.mods)[j])]
        
        if (!grepl("Z", full.mods[[j]][pars.tmp[1]])) {
          totnz <- totnz + res$freq[which(res$id == names(full.mods)[j])]
        }
        
        if (grepl("Z", full.mods[[j]][pars.tmp[1]])) {
          totz <- totz + res$freq[which(res$id == names(full.mods)[j])]
        }
      }

    }
      row <- which(rownames(pairwise) == pairs[ ,i][1])
      col <- which(colnames(pairwise) == pairs[ ,i][2])
      pairwise.f[row, col] <- tot
      pairwise.p[row, col] <- tot / length(mods)
      pairwise.f.nz[row, col] <- totnz
      pairwise.p.nz[row, col] <- totnz / length(mods)
      pairwise.f.jz[row, col] <- totz
      pairwise.p.jz[row, col] <- totz / length(mods)
  }
  
  pw.freq <- list(data.frame(pairwise.f), data.frame(pairwise.f.nz), data.frame(pairwise.f.jz))
  names(pw.freq) <- c("pairwise.all", "pairwise.nozero", "pairwise.justzero")
  pw.prob <- list(data.frame(pairwise.p), data.frame(pairwise.p.nz), data.frame(pairwise.p.jz))
  names(pw.prob) <- c("pairwise.all", "pairwise.nozero", "pairwise.justzero")
  
  ret <- list(res, zerotable, pw.freq, pw.prob)
  names(ret) <- c("models", "zero.params", "parameters.frequency", "parameters.probability")
  
  return(ret)
}





