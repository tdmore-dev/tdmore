#' TDMore
#'
#' Estimate individual parameters for population pharmacometrics models, and use them to make predictions and/or dose recommendations.
#'
#' @name tdmore
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tibble tibble
#' @docType package
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
utils::globalVariables(c("."))
.onAttach <- function(libname, pkgname) {
  if(interactive()) system.file("DISCLAIMER", package="tdmore") %>% readLines() %>% paste(collapse="\n") %>% packageStartupMessage() #nocov
}

utils::globalVariables(c("modName", "Rfile"))
