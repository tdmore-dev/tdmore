#' TDMore
#'
#' Estimate individual parameters for population pharmacometrics models, and use them to make predictions and/or dose recommendations.
#'
#' @name tdmore
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @docType package
NULL

.onAttach <- function(libname, pkgname) {
  if(interactive()) system.file("DISCLAIMER", package="tdmore") %>% readLines() %>% paste(collapse="\n") %>% packageStartupMessage() #nocov
}
