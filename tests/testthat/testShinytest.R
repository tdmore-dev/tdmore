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

tmp_lib <- tdmore::ensurePackagePresent("tdmore")

testpath <- testthat::test_path("../shinytest/")
message("Searching for shiny tests in ", testpath)
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
