#' @export
print.recommendation <- function(x, ...) {
  cat("Recommendation: \n")
  cat("A dose of `", x$dose, "` will hit the requested target of\n")
  print(x$target, row.names = F)
}

#' @export
summary.recommendation <- function(object, ...) {
  list(dose = object$dose,
       regimen = object$regimen)
}

#' @export
coef.recommendation <- function(object, ...) {object$dose}
