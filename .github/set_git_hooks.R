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
  text <- readLines(x)
  paste(text, collapse="\n")
}
script <- readFile(".github/pre-commit.R")
usethis::use_git_hook(hook = "pre-commit", script)
