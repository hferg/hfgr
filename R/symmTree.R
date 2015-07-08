#' symmTree
#'
#' Make a totally symmetrical tree of n tips with grafen branch lengths.
#' @param n The number of tips (must be apower of 4)
#' @param bls Branch lengths - either "grafen" or "equal". Defaults to "equal".
#' @export

symmTree <- function(n, bls = "equal") {

  # Make vector of tip labels.
  tips <- paste0("t", c(1:n))
  
  # Construct newick string - first make pairs of each of the tips.
  string <- tips

  while (length(string) >= 2) {
    tmp <- vector()

    for (i in seq(2, length(string), 2)) {
      tmp <- c(tmp, paste0("(", string[i - 1], ",", string[i], ")"))
    }
    
    string <- tmp
  }
  
  string <- paste0(string, ";")
  
  if (bls == "grafen") {
    tree <- compute.brlen(read.tree(text = string))
  } else if (bls == "equal") {
    tree <- compute.brlen(read.tree(text = string))
    tree$edge.length <- rep(1, nrow(tree$edge))
  }

  return(tree)
}

