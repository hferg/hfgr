#'
#' getPL
#' Gets the path length from the root to the position of a particular node.
#' @export
#' @name getPL

getPL <- function(tree, startnode = NA, node) {

  nodepath <- function(phy) .Call (seq_root2tip, phy$edge, Ntip(phy), phy$Nnode)

  # Get all paths for this tree
  allpaths <- nodepath(tree)
  # get all paths containing this node
  if(node > Ntip(tree)) {
    paths <- allpaths[grepl(node, allpaths)]
  } else {
    paths <- allpaths[[node]]
  }
  # If the terminal node we care about is not a terminal, we want to remove the branches after our node
  if(node > Ntip(tree) ) {
    paths <- lapply(paths, function(x) x <- x[x > Ntip(tree)])
    path <- unlist(unique(lapply(paths, function(x) x <- x[x <= node]))) 
  } else {
    path <- unlist(unique(paths))
  }
  # If the startnode is not NA, then we want to trim the paths. 
  if(!is.na(startnode)) {              
    while(path[1] != startnode) path <- path[-1]
    path <- path[-1] 
  } # Remove the start node [ otherwise it puts the leading branch in there too ]
  # get the path length
  branches <- which(tree$edge[,2] %in% path)
  distance <- sum(tree$edge.length[branches])
  return(distance)
}
