library(tdmore)
library(nlmixr)
library(ggplot2)
library(testthat)

context("Fully test the EBE routines")

# Load the default tdmore
m1 <- default_model %>% tdmore()
regimen <- data.frame(
  TIME=seq(0, 1)*24,
  AMT=5 #5mg
)
test_that("Crazy parameter estimates give fixed log likelihood", {
  res <- coef(m1)
  res[1] <- 1E6
  res[2] <- -1E6

  expect_error(
    estimate(m1,
             par = res, #crazy starting estimates
             observed = data.frame(TIME = 0:10, CONC = 10),
             regimen = regimen,
             covariates = c(WT = 70) )
  )

  crazyEstimate <- tdmorefit(m1, res=res, observed=data.frame(TIME=1,CONC=10),regimen=regimen,covariates=c(WT=70))
  expect_equal(logLik(crazyEstimate), as.numeric(NA))
})

test_that("Test all other ebe methods", {
  observed <- data.frame(TIME=10, CONC=0.04)
  fit <- estimate(m1, observed=observed,
                  regimen=regimen, covariates=c(WT=70))
  expect_known_output( print(fit), "tests/ref/print-tdmorefit-example" )
  expect_known_output( summary(fit), "tests/ref/summary-tdmorefit-example" )
  expect_known_value( confint(fit), "tests/ref/confint-tdmorefit-example" )
  expect_equal( model.frame(fit), observed )
  expect_known_value( fitted(fit), "tests/ref/fitted-tdmorefit-example")
  expect_equal( logLik(fit), 3.98965892617313 )
  a <- logLik(fit, type="pop" )
  b <- logLik(fit, type="pred" )
  expect_equal( logLik(fit), a+b )
  expect_equal( a, 0.682675762497937 )
  expect_equal( b, 3.30698316367519 )
  expect_error(logLik(fit, type="babla" ) )
  expect_true( is.tdmorefit(fit) )
  expect_true( !is.tdmorefit(m1) )
  # test predict
  expect_equal( predict(fit), data.frame(TIME=10, CONC=0.03142424) , tolerance=1E-6)
  # test predict with fixed values
  par <- coef(fit)
  expect_equal( predict(fit, parameters=par[1]*3 ), data.frame(TIME=10, CONC=0.03550358), tolerance=1E-6)
  # test predict with fixed values and se.fit
  predict(fit, se.fit=T)
  predict(fit, se.fit=T, parameters=par[1]*3)
  predict(fit, se.fit=T, level=NA)
  foo <- predict(fit, se.fit=T, level=NA, parameters=par[1]*3)
  expect_equal(foo$ECL, rep(-0.5314563, 100), tolerance=1E-7 )
})
