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

test_that("Test all other ebe methods", {
  observed <- data.frame(TIME=10, CONC=0.04)
  fit <- estimate(m1, observed=observed,
                  regimen=regimen, covariates=c(WT=70))
  expect_known_output( print(fit), "print-tdmorefit-example.txt" )
  expect_known_output( summary(fit), "summary-tdmorefit-example.txt" )
  expect_known_output( print(confint(fit)), "confint-tdmorefit-example.txt" )
  expect_equal( model.frame(fit), observed )
  expect_known_output( print(fitted(fit)), "fitted-tdmorefit-example.txt")
  expect_equal( logLik(fit), 3.98965892617313 )
  a <- logLik(fit, type="pop" )
  b <- logLik(fit, type="pred" )
  expect_equal( logLik(fit), a+b )
  expect_equal( a, 0.682675762497937 )
  expect_equal( b, 3.30698316367519 )
  expect_error( logLik(fit, type="babla" ) )
  expect_true( is.tdmorefit(fit) )
  expect_true( !is.tdmorefit(m1) )

  #How would a patient with higher bodyweight look? Can we re-estimate?
  fit2 <- estimate(fit, covariates=c(WT=75))
  expect_known_output( summary(fit), "summary-tdmorefit2-example.txt" )

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

test_that("We can use multiple methods and pick the best one", {
  observed <- data.frame(TIME=10, CONC=0.04)
  fit <- estimate(m1, observed=observed, regimen=regimen, covariates=c(WT=70))
  # Try with an invalid method...
  expect_error({
    expect_warning({
      estimate(fit, method="SANN")
    }, regexp="optimr: optim\\(\\) with bounds ONLY uses L-BFGS-B")
  }, )

  cg <- estimate(fit, lower=NULL, upper=NULL, method="CG")
  expect_equal( coef(cg), coef(fit), tolerance=1e-6 )

  assign("deps", 1e-4, optextras::optsp) #fix for grfwd in optextras
  rcgmin <- estimate(fit, lower=NULL, upper=NULL, method="Rcgmin")
  nlminb <- estimate(fit, lower=NULL, upper=NULL, method="nlminb")

  #SANN returns all 0, because it does not work
  #And it takes a long time to boot...
  #sann <- estimate(fit, lower=NULL, upper=NULL, method="SANN")
  #expect_true( all(coef(sann) == 0) )

  multi <- estimate(fit, lower=NULL, upper=NULL, method=c("nlminb", "Rcgmin"))
  expect_equal( coef(multi), coef(rcgmin) )
  expect_true( logLik(rcgmin) > logLik(nlminb) )

  #allmeth <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", #"nlm", #nlm has errors with Grad
  #             "nlminb",
  #             "Rcgmin", "Rvmmin", "hjn")


  multistart <- estimate(fit, multistart=TRUE)

  multistart <- estimate(fit, multistart=c(ECL=2, EV1=2))
  multistart <- estimate(fit, multistart=c(ECL=2))

  fitNoSe <- estimate(fit, se.fit=F)
  expect_true( all( diag(vcov(fitNoSe)) < 1e-8 ) )

  fixed <- estimate(fit, fix=c(EV1=1.2))
  expect_equal( coef(fixed)['EV1'] %>% unname , 1.2 )

  multiFix <- estimate(fit, fix=c(EV1=1.2), multistart=c(ECL=2))
  expect_equal(coef(multiFix), coef(fixed), tolerance=1e-5)
  expect_equal(logLik(multiFix), logLik(fixed))

  expect_error({
    estimate(fit, fix=c(EV1=1, EV1=2))
  }) # EV1 appears twice!
  expect_error({
    estimate(fit, fix=c(EV1=1, ECL=2))
  }) #all parameters fixed...
})
