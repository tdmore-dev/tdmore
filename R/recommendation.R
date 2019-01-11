#' Print a recommendation object.
#'
#' @param x a recommendation object
#' @param ... ignored
#'
#' @export
print.recommendation <- function(x, ...) {
  cat("Recommendation: \n")
  cat("A dose of `", x$dose, "` will hit the requested target of\n")
  print(x$target, row.names = F)
}

#' Summarize a recommendation object.
#'
#' @param object a recommendation object
#' @param ... ignored
#'
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
