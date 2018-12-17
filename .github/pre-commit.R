#!/c/Program Files/R/R-3.5.1/bin/x64/Rscript
##
## Script name: pre-commit.R
##
## Purpose of script:
## Verify good state of package before committing
##
## Author: Ruben Faelens
##
## Date Created: Wed Nov 14 10:24:13 2018
##
## Copyright (c) Ruben Faelens, 2018
## Email: ruben.faelens@gmail.com
##
## ---------------------------
##
## Notes:
##  The idea behind this hook is to never commit a *part* of your changes
##  and forget to commit something that should have been updated
##  We therefore have the following workflow:
##    what are you about to commit?
##    does this result in changes to dependent files?
##    do the file modification dates agree with dependencies?
##    ( are you committing these dependent files? ) -> no, as we assume this is a conscious decision
##
## This system will get confused if you made multiple modifications, but are only committing a part.
## E.g.
##   ChangeA: typo in source code file #1
##   ChangeB: typo in documentation for file #2, should also commit updated pkgdown site
##   The pre-commit hook for ChangeA will complain that the pkgdown site should be recompiled.
## ---------------------------
suppressMessages(
  library(tidyverse)
)

stat <- git2r::status()$staged %>% tibble(file=., change=names(.)) %>%
  unnest(file)

files <- list.files(recursive=TRUE) %>%
  tibble(file=.) %>%
  mutate(info = map(file, file.info) ) %>%
  {suppressWarnings(unnest(., info))} %>%
  mutate(staged=if(nrow(stat)==0) FALSE else file %in% stat$file )

## Function to convert a series of boolean expressions into a categorical number
categorize <- function(...) {
  booleans <- list(...)
  categories <- rep(NA, length(booleans[[1]]) )
  for(i in seq_along(booleans)) {
    categories[
      is.na(categories) &
        booleans[[i]]
      ] <- i
  }
  x <- factor(categories, levels=seq_along(booleans), labels=c(names(booleans)))
  as.character(x)
}

grepl <- function(x, pattern) {base::grepl(pattern=pattern, x=x)}
files <- files %>% mutate(category = categorize(
  source=startsWith(file, "R/"),
  tests=startsWith(file, "tests/"),
  instExamples=startsWith(file, "inst/examples"),
  vignettes=grepl(file, "^vignettes/.*Rmd$"),
  #note that built vignettes should not be in github, https://github.com/r-lib/devtools/issues/584
  #We simply use this as a shortcut for better dependency detection
  vignettesOut=grepl(file, "^vignettes/.*$"),
  book=startsWith(file, "bookdown/"),
  metadata=(file %in% c("DESCRIPTION","LICENSE")),
  bookResult=startsWith(file, "docs/book/"),
  pkgdown=(file == "_pkgdown.yml"),
  pkgdownResult= startsWith(file, "docs/"),
  man=startsWith(file, "man/"),
  NAMESPACE=file=="NAMESPACE",
#  READMERmd=(file=="README.Rmd"),
  READMEmd=(file == "README.md"), # this is used in github
  index.Rmd=(file=="index.Rmd"), # this file gets used for pkgdown::build_site
  Contributing=(file == "CONTRIBUTING.md"),
  Rproj=file=="tdmore.Rproj"
))



deps <- "category dep
NAMESPACE source
man source
bookResult book
bookResult vignettes
bookResult source
pkgdownResult pkgdown
pkgdownResult man
pkgdownResult vignettes
pkgdownResult source
pkgdownResult READMERmd
pkgdownResult metadata
vignettesOut vignettes
"
#Removed: (READMEmd is not generated anymore, but static developer information)
#READMEmd READMERmd
#READMEmd source
#
depsTable <- read.table(text=deps, header=TRUE)

categoryMtimes <- files %>%
  group_by(category) %>%
  arrange(mtime) %>% #from earliest to latest
  filter(row_number() == n()) %>% # pick the file with maximum mtime
  #summarize(maxMtime=max(mtime)) %>%
  merge(depsTable, all = TRUE) %>%
  filter(!is.na(mtime))
comparison <- categoryMtimes %>%
  filter(!is.na(dep)) %>%
  merge(categoryMtimes,
        by.x="dep", by.y="category",
        suffixes = c(".category", ".dep")) %>%
  ## The rule is as follows:
  ## If something changed, all the dependencies should have an mtime that is later
  ## Or otherwise put: the earliest mtime of something should always be LATER than the maximum mtime of dependencies
  ## The maxMTime of the category shows when it was last updated
  ## Using MinMTime does not make sense; e.g. the introduction of the book always stays the same
  mutate(OK = mtime.category >= mtime.dep)

problems <- comparison %>% filter(!OK) %>%
  select(dep, file.dep, mtime.dep, category, file.category, mtime.category)

## TODO: Check the currently staged files. If they are not contained in
## this problem list (as source or dependent files), then we could still commit...

if(nrow(problems) > 0) {
  # TODO: add some nicer output by using the crayon package
  # https://github.com/r-lib/crayon#usage
  message("Files <", paste(problems$category, collapse=", "), "> were not updated when their dependencies <", paste(unique(problems$dep), collapse=", "), "> changed.")
  problems %>%
    merge(files, by.x="dep", by.y="category") %>%
    filter(mtime.dep > mtime.category) %>%
    select(category, mtime.category, file.category, dep, mtime, file) %>%
    print

  message("Executing the required update functions... ")
  if(any(c("man", "NAMESPACE") %in% problems$category)) {
    if("man" %in% problems$category)
      do.call(file.remove, list(list.files("man/", full.names = TRUE)))
    devtools::document()
  }
  if("READMEmd" %in% problems$category) devtools::build_readme(quiet=TRUE)
  if("pkgdownResult" %in% problems$category) {
    message("Building site...")
    devtools::build_site(quiet=TRUE)
  }
  if("bookResult" %in% problems$category) source("bookdown/create_bookdown.R")

  message("Review the modified files, stage them, and commit again")
  stop(call. = FALSE)
}

