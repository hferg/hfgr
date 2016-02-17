#' simuDatEb
#'
#' Simulate data over a phylogeny, or sub section of a phylogeny, in order to reflect variable rates.
#' Based on adding a constant to the tip data branch by branch and weighting by branch length and does
#' not need to stretch the phylogeny. This is good for exponential rate decay, but pretty shoddy for delta
#' at the moment - what it CAN do is do a good approximation of the delta relationship, but t isn't functioning
#' as half life at the moment.

#' @param tree The tree that data is simulated on.
#' @param node The node at which to make the changes.
#' @param a Initial rate increase before decay. If using the power decay this is, essentially, delta.
#' @param t In "pwr" (delta) mode this is ROUGHLY the delta parameter. In "exp" earyl burst mode, it is the number of half lives for the decay to take place over. Essentially a high number means the rate decays rapidly, and a low one means it decays slowly. This is used rather than a rate parameter, since decay is relative to branch lengths, but half lives are not.
#' @param sig The standard deviation of the random component of the BM simulation. NOT sigma squared.
#' @param decay Either "exp" (exponential) or "pwr" (power). The shape of the decay from the initial rate.
#' @param dat If you want to use some data you already have and change it this way, put it here.
#' @export

simuDatEB <- function(tree, node, sig, a, t, decay = "exp", dat = NULL) {
  
  if (isDefined(dat)) {
    dat <- dat
  } else {
    dat <- rTraitCont(tree, model = "BM", sigma = sig)
  }
  
  transdat <- dat  

  descs <- c(getDescs(tree, node), node)
  internal <- descs[descs > length(tree$tip.label)]
  tips <- descs[descs < length(tree$tip.label)]
  edges <- tree$edge[tree$edge[ ,2] %in% c(internal, tips), ]

  heights <- nodeHeights(tree)
  heights <- heights[tree$edge[ ,2] %in% descs, ]
  age <- max(heights)

  transnode <- nodeheight(tree, node)
  transtime <- age - transnode



  if (decay == "exp") {
    hl <- transtime / t
    r <- log(2) / hl
    fun <- function(x) a * (exp(1) ^ (-r * x))
  } else if (decay == "pwr") {
    if (t >= 1) {
      fun <- function(x) x ^ (1 / t)
    } else if (t < 1) {
      fun <- function(x) x ^ -(1 / t)
    }

  }
  
  for (i in 1:nrow(edges)) {
    tps <- getDescs(tree, edges[i, 2])
    tps <- tps[tps < length(tree$tip.label)]
    num <- sample(2, 1)
    
    start <- transtime - (age - nodeheight(tree, edges[i, 1]))
    end <- transtime - (age - nodeheight(tree, edges[i, 2]))
    
    m <- integrate(fun, start, end)$value
    
    if (num == 1) {
      transdat[tps] <- transdat[tps] + (m * sig)
    } else {
      transdat[tps] <- transdat[tps] - (m * sig)
    }
    
  }

  res <- list(original = dat, transformed = transdat)
  return(res)
}

