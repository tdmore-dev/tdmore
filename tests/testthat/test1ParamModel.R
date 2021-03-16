library(tdmore)
library(RxODE)
library(nlmixr)
library(ggplot2)
library(testthat)

set.seed(0)

# Load an example model
m1 <- getModel("example_1param")

# Predicting new data
regimen <- data.frame(
  TIME=c(0, 8, 16),            # Every 8 hour and for 1 day, an injection is given
  AMT=c(1000, 1000, 1000),     # 1g is administered
  RATE=c(1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
)

# Estimate and plot IPRED
observed <- data.frame(TIME=c(9, 16), CONC=c(30, 8))
ipred <- m1 %>% tdmore:::estimate(observed = observed, regimen = regimen)
test_that("IRES and IWRES make sense", {
  expect_snapshot_output(residuals(ipred, weighted=FALSE))
  expect_snapshot_output(residuals(ipred, weighted=TRUE))
})


test_that("Prediction results makes sense", {
  z1 <- plot(ipred)
  expect_doppelganger("ipred_1cmt", z1)
  z2 <- tdmore:::autoplot.tdmoreprofile(profile(ipred))
  expect_doppelganger("ipred_1cmt_profile", z2)
  expect_equal(round(ipred$res, digits=3), c(EV1=0.578))
})
