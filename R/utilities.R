#' Convert a numeric vector/list into a diagonal matrix.
#' If the vector/list is named, names will be used as the column and row names of the resulting diagonal matrix.
#'
#' @param vector a numeric vector/list
#'
#' @return a diagonal matrix
#' @export
vectorToDiagonalMatrix <- function(vector) {
  diagMatrix <- diag(x = vector)
  if (!is.null(names(vector))) {
    colnames(diagMatrix) <- names(vector)
    rownames(diagMatrix) <- names(vector)
  }
  return(diagMatrix)
}

#' Melt prediction results.
#'
#' @param x dataframe returned by the predict function
#'
#' @return the melted dataframe
meltPredictions <- function(x, se=FALSE) {
  measure.vars <- colnames(x)
  measure.vars <- measure.vars[measure.vars != "TIME"]
  vars <- measure.vars
  if(se) {
    for(i in c(".upper", ".lower")) vars <-c(vars, paste0(measure.vars, i))
  }
  tmp <- reshape::melt(x, id.vars="TIME")
  if(se) {
    result <- tmp[ tmp$variable %in% measure.vars , ]
    result$value.upper <- tmp$value[ tmp$variable %in% paste0(measure.vars, ".upper") ]
    result$value.lower <- tmp$value[ tmp$variable %in% paste0(measure.vars, ".lower") ]
  } else {
    result <- tmp
  }
  return(result)
}

#' Compute the TMax value based on the specified regimen and observed data.
#'
#' @param regimen the specified regimen
#' @param observed the observed dataframe
#'
#' @return the tmax value, numeric
computeTmax <- function(regimen, observed=NULL) {
  values <- c(0)
  if (!is.null(regimen)) {
    values <- c(values, regimen$TIME, regimen$TIME + regimen$ADDL * regimen$II)
  }
  if (!is.null(observed)) {
    values <- c(values, observed$TIME)
  }
  return(max(values, na.rm=TRUE))
}
