#' smartBind 
#'
#' A function that will rbind vectors of different lengths and return a matrix, provided each vector element is named.
#' @param A bunch of vectors.
#' @export
#' @examples
#' do.call(smartRbind, list.of.vectors)

smartBind <- function (...) {
  # from GSee http://stackoverflow.com/questions/17308551/do-callrbind-list-for-uneven-number-of-column
  dargs <- list(...)

  if (!all(vapply(dargs, is.vector, TRUE))) 
      stop("all inputs must be vectors")

  if (!all(vapply(dargs, function(x) !is.null(names(x)), TRUE))) 
      stop("all input vectors must be named.")

  all.names <- unique(names(unlist(dargs)))
  out <- do.call(rbind, lapply(dargs, `[`, all.names))
  colnames(out) <- all.names
  out
}


