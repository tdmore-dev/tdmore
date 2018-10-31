#' Instantiate a new error model.
#'
#' @param var the related variable name
#' @param add additive residual error, as stdev
#' @param prop proportional residual error, as stdev
#' @param exp exponential residual error, as stdev. The exponential error cannot be used in conjunction with the additive or proportional error
#'
#' @export
errorModel <- function(var, add=0, prop=0, exp=0) {
  assertthat::is.string(var)

  return(structure(
    list(
      var = var,
      add = add,
      prop = prop,
      exp = exp
    ),
    class = c("error_model")
  ))
}
