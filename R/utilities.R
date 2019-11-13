#' Convert a numeric vector/list into a diagonal matrix.
#' If the vector/list is named, names will be used as the column and row names of the resulting diagonal matrix.
#'
#' @param vector a numeric vector/list
#'
#' @return a diagonal matrix
vectorToDiagonalMatrix <- function(vector) {
  if(length(vector)==1) {
    diagMatrix <- matrix(data=vector, nrow=1, ncol=1)
  } else {
    diagMatrix <- diag(x = vector)
  }
  if (!is.null(names(vector))) {
    colnames(diagMatrix) <- names(vector)
    rownames(diagMatrix) <- names(vector)
  }
  return(diagMatrix)
}

#' Melt prediction results.
#'
#' @param x dataframe returned by the predict function
#' @param se if FALSE, melts the data.frame as expected
#' if TRUE, assumes every output has 4 columns in the data.frame: XXX1, XXX1.median, XXX1.lower, XXX1.upper
#' And melts these into a variable (with values XXX1, XXX2, ...) and 4 columns: value, value.median, value.lower and value.upper
#'
#' @return the melted dataframe
#' @keywords internal
meltPredictions <- function(x, se=FALSE) {
  gathered <- x %>%
    tidyr::gather(key = "variable", value="value", -.data$TIME)
  if(!se) return(gathered)

  gathered %>% tidyr::extract(.data$variable, into=c("variable", "suffix"),
                              regex="^(.*?)(|\\.lower|\\.upper|\\.median)$") %>%
    dplyr::mutate(suffix=paste0("value", .data$suffix) ) %>%
    tidyr::spread(key=.data$suffix, value=.data$value)
}
