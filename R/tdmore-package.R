#' TDMore
#'
#' Estimate individual parameters for population pharmacometrics models, and use them to make predictions and/or dose recommendations.
#'
#' @name tdmore
#' @importFrom magrittr %>%
#' @docType package
NULL

.onAttach <- function(libname, pkgname) {
  if(interactive()) packageStartupMessage("Welcome to TDMore, more info on https://tdmore-dev.github.io/tdmore.") #nocov
}
