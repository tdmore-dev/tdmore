library(tdmore)
library(RxODE)
library(nlmixr)
library(ggplot2)
library(testthat)

context("Test a 1-parameter-model")

set.seed(0)

# Load the Meropenem model
tdmore <- getModel("meropenem_1param")

# Predicting new data
regimen <- data.frame(
  TIME=c(0, 8, 16),            # Every 8 hour and for 1 day, an injection is given
  AMT=c(1000, 1000, 1000),     # 1g is administered
  RATE=c(1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
)

# Estimate and plot IPRED
observed <- data.frame(TIME=c(9, 16), CONC=c(30, 8))
ipred <- tdmore %>% estimate(observed = observed, regimen = regimen)
test_that("IRES and IWRES make sense", {
  expect_known_value(residuals(ipred, weighted=FALSE), "1-param-IRES")
  expect_known_value(residuals(ipred, weighted=TRUE), "1-param-IWRES")
})

plot(ipred)
test_that("Prediction results makes sense", {
  expect_doppelganger("ipred_1cmt", plot(ipred))
  expect_doppelganger("ipred_1cmt_profile", plot(profile(ipred)))
  expect_equal(round(ipred$res, digits=3), c(EV1=0.578))
})
