#' asymmTree
#'
#' Simulates a completely unbalanced tree with n tips, and branch lengths calculated according to Grafen.
#' @param n The number of tips for the tree to have.
#' @export

asymmTree <- function(n) {

  # Make vector of tip labels.
  tips <- paste0("t", c(1:n))
  
  # Construct newick string.

  string <- paste0("(", tips[1], ",", tips[2], ")")
  
  tips2 <- tips[3:length(tips)]
  
  for (i in tips2) {
    string <- paste0(string, ",", i, ")")
  }
  
  # Add open brackets to the front to close all brackets (count ")" and subtract 1), and semicolon to end.
  cb <- gsub(")", "", string)
  reps <- (nchar(string) - nchar(cb)) - 1
  
  string <- paste0(paste(rep("(", reps), collapse = ""), string, ";")
  
  tree <- compute.brlen(read.tree(text = string))
  
  return(tree)
}

