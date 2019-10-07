#' TDMore
#'
#' Estimate individual parameters for population pharmacometrics models, and use them to make predictions and/or dose recommendations.
#'
#' @name tdmore
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @docType package
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
utils::globalVariables(c("."))

loadModel <- function(modName, Rfile) {
  if(interactive()) cat("Loading model ", modName, " ...\n")
  result <- source(Rfile, keep.source = TRUE)
  result$value
}

.onAttach <- function(libname, pkgname) {
  if(interactive()) system.file("DISCLAIMER", package="tdmore") %>% readLines() %>% paste(collapse="\n") %>% packageStartupMessage() #nocov
}

utils::globalVariables(c("modName", "Rfile"))

## Lazy load all models
.onLoad <- function(libname, pkgname) {
  modelDir <- system.file("models", package=pkgname)
  files <- dir(modelDir)
  for(i in tools::file_path_sans_ext(files)) {
    e1 <- new.env()
    e1$Rfile = file.path(modelDir, paste0(i, ".R") )
    e1$modName = i
    delayedAssign(x=i,
                  value={
                    loadModel(modName, Rfile)
                  },
                  eval.env=e1,
                  assign.env=parent.env(environment()))
  }
}
