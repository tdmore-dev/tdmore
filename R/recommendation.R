#' Print a recommendation object.
#'
#' @param object a recommendation object
#'
#' @return print the recommended dose
#' @export
print.recommendation <- function(object, ...) {
  print(object$dose)
}

#' Summarize a recommendation object.
#'
#' @param object a recommendation object
#'
#' @return print the recommended dose
#' @export
summary.recommendation <- function(object, ...) {
  list(dose = object$dose,
       regimen = object$regimen)
}

#' The obtained parameter values for a recommendation.
#'
#' @param object a recommendation object
#' @param ... ignored
#'
#' @return a named numeric vector
#' @export
coef.recommendation <- function(object, ...) {object$dose}
