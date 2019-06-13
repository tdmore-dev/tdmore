library(shinytest)
library(testthat)

context("Test shinyProfile")

test_that("shinyProfile behaves as expected", {
  skip_on_cran()

  appDir <- "../testShinyProfile/"
  expect_true(file.exists(file.path(appDir, "app.R")))

  ## Some caveats
  ## 1) PhantomJS uses the QT library to render
  ## Different platforms will render differently (e.g. buttons are rendered differently on Linux vs Windows)
  ## 2) PhantomJS will compare the content of plots
  ## Because of #1, the plotting area may have a different size.
  ## Also, fonts may be different as well, which results in a different SHA1 hash for the PNG image.
  ## We recommend to *ignore* all plots, but instead to do one of the following
  ##    a) export the source data of the plot using exportTestValues()
  ##    b) use a htmlwidget (like plotly) to render the plot

  expect_pass(
    testApp(appDir=appDir, compareImages=interactive(), quiet=interactive())
  )
})
