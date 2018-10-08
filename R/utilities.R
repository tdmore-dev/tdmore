#' Convert a numeric vector/list into a diagonal matrix.
#' If the vector/list is named, names will be used as the column and row names of the resulting diagonal matrix.
#'
#' @param vector a numeric vector/list
#'
#' @return a diagonal matrix
#' @export
convertVectorToDiag <- function(vector) {
  diagMatrix <- diag(x = vector)
  if (!is.null(names(vector))) {
    colnames(diagMatrix) <- names(vector)
    rownames(diagMatrix) <- names(vector)
  }
  return(diagMatrix)
}
