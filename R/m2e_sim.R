#' simuleteFlips
#' 
#' Simulates the outcome of opposed duels.
#' @param attack Attackers attack value
#' @param defense Defenders resist value
#' @param TN Target number (for simple duel - ignores any defense stuff)
#' @param suit Is a specific suit needed?
#' @param aplus Number of positive twists for attacker
#' @param aminus Number of negative twists for attacker
#' @param dplus Number of positive twists for defender
#' @param dminus Number of negative twists for defender
#' @param dtrack Currently ignored
#' @name simulateFlips
#' @export

simulateFlips <- function(attack, defense, TN = NULL, suit = FALSE, aplus = 0, aminus = 0, dplus = 0, dminus = 0, dtrack = NULL) {
  # TO_THINK
  #   Should this include some attempt to make the deck smaller based on what has come out of it?
  # Build the deck, include both jokers.
  require(ggplot2)

  deck <- c(rep(c(1:13), 4), 0, 14)
  nsim <- 10000
  attackerfate <- aplus - aminus
  defenderfate <- dplus - dminus
  
  attackerresults <- vector(mode = "numeric", length = nsim)
  defenderresults <- vector(mode = "numeric", length = nsim)

  if (attackerfate < 0) {
    attackerflips <- abs(attackerfate) + 1
    print(paste("Attacker flips", attackerflips, "cards"))

    for (i in 1:nsim) {
      attackerresults[i] <- attack + min(sample(deck, attackerflips, replace = FALSE))
    }
  } else {
    attackerflips <- abs(attackerfate) + 1
    print(paste("Attacker flips", attackerflips, "cards"))    

    for (i in 1:nsim) {
      attackerresults[i] <- attack + max(sample(deck, attackerflips, replace = FALSE))
    }
  }

  if (!is.null(TN)) {
    defenderresults <- rep(TN, nsim)
  } else {

    if (defenderfate < 0) {
      defenderflips <- abs(defenderfate) + 1
      print(paste("Defender flips", defenderflips, "cards"))

      for (i in 1:nsim) {
        defenderresults[i] <- defense + min(sample(deck, defenderflips, replace = FALSE))
      }
    } else {
      # choose big
      defenderflips <- abs(defenderfate) + 1
      print(paste("Defender flips", defenderflips, "cards"))

      for (i in 1:nsim) {
        defenderresults[i] <- defense + max(sample(deck, defenderflips, replace = FALSE))
      }
    }
  }

  diffs <- attackerresults - defenderresults
  
  psuccess <- round(sum(diffs >= 0) / length(diffs), 2)


  # make a nice histogram to show the attacker success distribution, and the probability of
  # success, and then the average damage flip as well.
  plot.cols <- c("orangered", "dodgerblue")                   
  names(plot.cols) <- c("fail", "success")
  bwidth <- 3.5 * sd(diffs) * length(diffs)^-(1/3)
  p <- data.frame(p = diffs, z = NA)
  colnames(p) <- c("dffs", "z")
  p$z[which(diffs < 0)] <- "fail"
  p$z[which(diffs >= 0)] <- "success"

  density_plot <- ggplot() +
          geom_density(data = p, aes(x = dffs, fill = "z"), alpha = 0.3) +
          geom_vline(data = p, aes(xintercept = 0), colour = "blue", linetype = "dashed", size = 1)
  dpb <- ggplot_build(density_plot)
  x1w <- max(which(dpb$data[[1]]$x <= 0))
  x2w <- max(which(dpb$data[[1]]$x > 0))
  x1m <- max(which(dpb$data[[1]]$x <= 1))
  x2m <- max(which(dpb$data[[1]]$x > 1))
  
  if (any(dpb$data[[1]]$x > 6)) {
    x1h <- max(which(dpb$data[[1]]$x <= 6))
    x2h <- max(which(dpb$data[[1]]$x > 6))
  } else {
    x1h <- x2m
    x2h <- x2m
  }

  if (any(dpb$data[[1]]$x > 11)) {  
    x1e <- max(which(dpb$data[[1]]$x <= 11))
    x2e <- max(which(dpb$data[[1]]$x > 11))  
  } else {
    x1e <- x2h
    x2e <- x2h
  }

  ret <- density_plot +
          geom_area(data = data.frame(x = dpb$data[[1]]$x[x1w:x2w],
                                      y = dpb$data[[1]]$y[x1w:x2w]),
            aes(x = x, y = y), fill = "dodgerblue", alpha = 0.3) +
          geom_area(data = data.frame(x = dpb$data[[1]]$x[x1m:x2m],
                                      y = dpb$data[[1]]$y[x1m:x2m]),
            aes(x = x, y = y), fill = "dodgerblue", alpha = 0.3) +          
          geom_area(data = data.frame(x = dpb$data[[1]]$x[x1h:x2h],
                                      y = dpb$data[[1]]$y[x1h:x2h]),
            aes(x = x, y = y), fill = "dodgerblue", alpha = 0.3) +          
          geom_area(data = data.frame(x = dpb$data[[1]]$x[x1e:x2e],
                                      y = dpb$data[[1]]$y[x1e:x2e]),
            aes(x = x, y = y), fill = "dodgerblue", alpha = 1) +                                
          ggtitle(paste0("p of success = ", psuccess))

  return(ret)
}



22 - 31