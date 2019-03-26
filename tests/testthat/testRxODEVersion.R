context("Test RxODE version running the tests")

### Newer versions of RxODE give slightly different results from current versions
### This results in test failures on very minor points.

library(RxODE)
test_that("RxODE version matches the one used to record the tests", {
  expect_known_value(sessionInfo()$otherPkgs$RxODE$Version, "rxode-version.txt" )
  expect_known_value(RxODE::rxVersion(), "rxode-dparser-version.txt" )
})
