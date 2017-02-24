#' phyloHeights
#'
#' Essentially the same function as heights.phylo from geiger
#'
#' @name phyloHeights

phyloHeights <- function(tree){
  tree <- reorder(tree, "postorder")
  n <- length(tree$tip.label)
  n.node <- tree$Nnode
  xx <- numeric(n + n.node)
  for (i in nrow(tree$edge):1) { 
    xx[tree$edge[i, 2]] <- xx[tree$edge[i, 1]] + tree$edge.length[i] 
  }
  

  root <- ifelse(is.null(tree$root.edge), 0, tree$root.edge)
  labs <- c(tree$tip.label, tree$node.label)
  depth <- max(xx)
  tt <- depth - xx
  idx <- 1:length(tt)
  dd <- tree$edge.length[idx]
  mm <- match(1:length(tt), c(tree$edge[, 2], Ntip(tree) + 1))
  dd <- c(tree$edge.length, root)[mm]
  ss <- tt + dd
  res <- cbind(ss, tt)
  rownames(res) <- idx
  colnames(res) <- c("start", "end")
  res <- data.frame(res)
  return(res)
}
