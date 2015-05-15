#' is.defined
#'
#' Internal function that tests if a variable is defined (i.e. is not NULL).
#' @param x The variable or object you want to test the existence of.

isDefined <- function(x) !is.null(x)
