library(shinytest)
library(testthat)

context("Test shinyProfile")

test_that("shinyProfile behaves as expected", {
  skip_on_cran()

  appDir <- "../testShinyProfile/"
  expect_true(file.exists(file.path(appDir, "app.R")))

  expect_pass(
    testApp(appDir=appDir, compareImages=interactive(), quiet=interactive())
  )
})
