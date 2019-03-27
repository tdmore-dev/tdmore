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
library(vdiffr)

rxOdeVersion <- function() {
  paste0( "RxODE-",sessionInfo()$otherPkgs$RxODE$Version )
}

## Change the path for different versions of RxODE...
expect_doppelganger_RxODE <- function(title, fig, path=NULL, ...) {
  vdiffr::expect_doppelganger(title, fig,
                              path=paste0("../",rxOdeVersion()),
                              ...)
}

expect_known_value_RxODE <- function(object, file, ...) {
  file <- paste0(rxOdeVersion(), "/", file)
  testthat::expect_known_value(object, file, ...)
}
