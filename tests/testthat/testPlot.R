library(tdmore)
library(nlmixr)
library(ggplot2)
library(testthat)

set.seed(0)

# Load the default tdmore
tdmore <- getModel("default")

regimen <- data.frame(
  TIME=seq(0, 1)*24,
  AMT=5 #5mg
)
covariates = c(WT=70)

# Default tdmore plot
test_that("plot() with tdmore object produces typical value plot", {
  z1 <- plot(tdmore, regimen, covariates=covariates) +
    coord_cartesian(ylim=c(0, 0.1))
  expect_doppelganger("default-tdmore-plot", z1)
})

# Create the observed and covariates dataframe
observed <- data.frame(TIME=2, CONC=0.040)

# Compute PRED
pred <- tdmore %>% estimate(regimen = regimen, covariates=covariates)
test_that("estimation without observed values produces ETA=0", {
  expect_equal(pred$res, c(ECL=0.0, EV1=0.0))
})

test_that("Plotting population prediction", {
  z1 <- plot(pred, covariates=covariates, newdata=seq(0, 48, by=0.1)) +
    coord_cartesian(ylim=c(0, 0.12))
  expect_doppelganger("pred-tdmore-plot", z1)
})

# Compute IPRED
test_that("estimation with observed values produced individual estimate", {
  ipred <- tdmore %>% estimate(observed = observed, regimen = regimen, covariates = covariates)
  expect_equal(round(ipred$res, digits=4), c(ECL=0.0336, EV1=0.1175))

  p1 <- plot(ipred, newdata=data.frame(TIME=seq(0, 48, by=0.1), CONC=NA)) +
    coord_cartesian(ylim=c(0, 0.1))
  expect_doppelganger("ipred-tdmore-plot", p1)
})

observed <- data.frame(TIME=c(2,26), CONC=c(0.040, 0.0675))
covariates <- data.frame(TIME=0, WT=70)
ipred <- tdmore %>% estimate(observed = observed, regimen = regimen, covariates = covariates)
ipred$res <- signif(ipred$res, 6)
ipred$varcov <- signif(ipred$varcov, 4)

## Different systems will have different digits after 1e-6...
## expect_known_output(dput(ipred), "testPlot-ipred")

set.seed(0)
test_that("Find dose example and test", {
  p2 <- plot(ipred, newdata=data.frame(TIME=seq(0, 48, by=0.1), CONC=NA)) +
    coord_cartesian(ylim=c(0, 0.12))
  expect_doppelganger("ipred2-tdmore-plot", p2)
})

set.seed(0)
test_that("plot() with numeric vector as newdata", {
  p2 <- plot(ipred, newdata=seq(0, 48, by=0.1)) +
    coord_cartesian(ylim=c(0, 0.12))
  expect_doppelganger("ipred2-tdmore-plot2", p2)
})
