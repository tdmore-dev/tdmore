##
## Script name: set_git_hooks.R
##
## Purpose of script:
## Set up git pre-commit hooks to make life easier
##
## Author: Ruben Faelens
##
## Date Created: Wed Nov 14 10:07:56 2018
##
## Copyright (c) Ruben Faelens, 2018
## Email: ruben.faelens@gmail.com
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

readFile <- function(x) {
  RscriptCall <- callr:::setup_rscript_binary_and_args(list())
  RscriptCall <- paste("#!", RscriptCall$bin, paste="")
  text <- readLines(x)
  text <- paste(text, collapse="\n")
  text <- paste(c(RscriptCall, text), paste="", collapse="\n")
}
script <- readFile(".github/pre-commit.R")
usethis::use_git_hook(hook = "pre-commit", script)
