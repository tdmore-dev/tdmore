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
#' @param se if FALSE, melts the data.frame as expected
#' if TRUE, assumes every output has 4 columns in the data.frame: XXX1, XXX1.median, XXX1.lower, XXX1.upper
#' And melts these into a variable (with values XXX1, XXX2, ...) and 4 columns: value, value.median, value.lower and value.upper
#'
#' @return the melted dataframe
meltPredictions <- function(x, se=FALSE) {
  #tmp <- reshape::melt(x, id.vars="TIME")
  if(se) {
    vars <- colnames(x)
    vars <- vars[vars != "TIME"]

    trueVars <- paste0(vars, ".median") %in% vars &
      paste0(vars, ".lower") %in% vars &
      paste0(vars, ".upper") %in% vars
    trueVars <- vars[trueVars]

    tmp <- reshape::melt(x, id.vars="TIME")

    result <- tmp[ tmp$variable %in% trueVars, ]
    result$value.median <- tmp$value[ tmp$variable %in% paste0(trueVars, ".lower") ]
    result$value.lower <- tmp$value[ tmp$variable %in% paste0(trueVars, ".lower") ]
    result$value.upper <- tmp$value[ tmp$variable %in% paste0(trueVars, ".upper") ]
  } else {
    result <- reshape::melt(x, id.vars="TIME")
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
