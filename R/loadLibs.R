#' loadLibs
#' 
#' Just a shortcut function to load all of the usual libraries at once.
#' @export
#' @name loadLibs

loadLibs <- function() {
  library(ape)
  library(ggplot2)
  library(geiger)
  library(phytools)
  library(reshape2)
  library(BTRTools)
  library(caper)
  library(TreeSim)
  library(ggtree)
  library(grid)
  library(pbapply)
  library(paleotree)
}