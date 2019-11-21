## This file searches for all directories within "shinytest"
## and executes an appropriate testApp
##
## Some caveats
## 1) PhantomJS uses the QT library to render
## Different platforms will render differently (e.g. buttons are rendered differently on Linux vs Windows)
## 2) PhantomJS will compare the content of plots
## Because of #1, the plotting area may have a different size.
## Also, fonts may be different as well, which results in a different SHA1 hash for the PNG image.
## We recommend to *ignore* all plots, but instead to do one of the following
##    a) export the source data of the plot using exportTestValues()
##    b) use a htmlwidget (like plotly) to render the plot

library(shinytest)
library(testthat)

context("Perform all shinytest tests")



tmp_lib <- tempfile("R_LIBS")
dir.create(tmp_lib)
on.exit({
  unlink(tmp_lib, recursive = TRUE)
}, add=T)

## Is the tdmore library available on the search path? If yes, then use that one
## If not, then install the tdmore library in a separate path, because shiny is executed in a separate process
pkg <- "tdmore"
pkgPat <- path.package(pkg, quiet=T)
if(is.null(pkgPat) || pkg %in% devtools::dev_packages() || startsWith(pkgPat, Sys.getenv("R_LIBS_USER")) ) {
  #Either the package is not available
  #or it is currently in the dev packages
  #or we are about to load the (outdated) installed version!
  message("Installing temporary version of tdmore for use with shinytest")
  pkg <- devtools::as.package(".")
  utils::install.packages(repos = NULL,
                          lib = tmp_lib,
                          pkg$path,
                          type = "source",
                          INSTALL_opts = c("--example",
                                           "--install-tests",
                                           "--with-keep.source",
                                           "--with-keep.parse.data",
                                           "--no-multiarch"),
                          quiet = T)
}


onCi <- !isTRUE(as.logical(Sys.getenv("CI")))
testpath <- rprojroot::find_package_root_file("tests/shinytest")
tests <- dir(testpath)

for(appDir in tests) {
  appDir <- appDir
  test_that(paste0("shinytest for ", appDir, "..."), {
    expect_pass(
      withr::with_libpaths(tmp_lib, {
        appDir <- file.path(testpath, appDir)
        testApp(appDir=appDir, compareImages=!testthat::is_testing())
      }, action="prefix")
    )
  })
}
