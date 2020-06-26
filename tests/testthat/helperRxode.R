##
## Script name: helperRxode.R
##
## Purpose of script:
## RxODE produces different results between 0.8.0-9 and 0.8.1
## We offer special methods to store results when using the development version of RxODE
##
## Author: Ruben Faelens
##
## Date Created: Wed Mar 27 14:50:54 2019
##
## Copyright (c) Ruben Faelens, 2019
## Email: ruben.faelens@gmail.com
##
## ---------------------------
##
## Notes:
## notes
##
## ---------------------------
library(testthat)
library(RxODE)

myVersion <- NULL
rxOdeVersion <- function() {
  if(!is.null(myVersion)) return(myVersion)

  testRoot <- rprojroot::find_testthat_root_file("ref")
  available <- dir(testRoot, pattern="RxODE-*")
  candidate <- paste0( "RxODE-",sessionInfo()$otherPkgs$RxODE$Version )
  if(candidate %in% available) {
    myVersion <<- candidate
  } else {
    ## Use the results from the previous version
    all <- sort( c(available, candidate) )
    i <- which( all == candidate )
    if(i == 1) stop("RxODE version too old, no available testthat repository!")
    myVersion <<- all[ i-1 ] #previous version
    warning("Using older version for RxODE testing repository: current=", candidate, ", reference=", myVersion )
  }
  return(myVersion)
}

rxOdePath <- function(...) {
  testthat::test_path("ref", rxOdeVersion(), ...)
}

setup({
  message("Tests using RxODE version: ", rxOdeVersion())
  message("Reference path: ", normalizePath(rxOdePath()))
})

## Change the path for different versions of RxODE...
expect_doppelganger <- function(title, fig, path=rxOdePath()) {
  # expect_doppelganger changes the path as test_path("..", "figs", path)
  # too bad; vdiffr::manage_cases() also expects the files in that path
  #stop("THIS FUNCTION IS CALLED with title=", title, ", fig=(notShown), and path=", path)
  vdiffr::expect_doppelganger(title=title, fig=fig, path=path, verbose=TRUE)
}

expect_known_value <- function(object, file, ...) {
  testthat::expect_known_value(object, rxOdePath(file), ...)
}

expect_known_output <- function(object, file, ...) {
  testthat::expect_known_output(object, rxOdePath(file), ...)
}
