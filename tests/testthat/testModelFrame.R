library(testthat)

m1 <- getModel("example")
test_that("model.frame correctly returns values", {
  expect_error( model.frame(m1) )
  expect_equal(
    model.frame(m1, data=0:10),
    data.frame(TIME=0:10, CONC=as.numeric(NA))
  )

  expect_equal(
    model.frame(m1, data=data.frame(TIME=0:10, CONC=13)),
    data.frame(TIME=0:10, CONC=13)
  )

})

test_that("model.frame correctly returns SE values", {
  expect_error( model.frame(m1, se=TRUE) )
  expect_equal(
    model.frame(m1, data=0:10, se=TRUE),
    data.frame(TIME=0:10, CONC=as.numeric(NA), CONC.lower=as.numeric(NA), CONC.upper=as.numeric(NA))
  )

  expect_equal(
    model.frame(m1, data=data.frame(TIME=10, CONC=20), se=TRUE),
    data.frame(TIME=10, CONC=20, CONC.lower=5.4570672347128, CONC.upper=34.5429327652872)
  )

  set.seed(1234) #set seed, to ensure we do not encounter an alpha error
  N <- 2000
  CONC <- replicate(N, model.frame(m1, data=data.frame(TIME=10, CONC=20), se=TRUE, level=NA)$CONC)
  sd <- m1$res_var[[1]]$prop
  sem <- sd / sqrt(N)
  sesd <- sd / sqrt(N-1) #standard error of sd, see https://en.wikipedia.org/wiki/Standard_deviation#Sample_standard_deviation_of_metabolic_rate_of_northern_fulmars

  expect_true( dplyr::between( mean(CONC/20)-1, qnorm(0.05, 0, sem), qnorm(0.95, 0, sem) ) )
  expect_true( dplyr::between( sd(CONC/20), qnorm(0.05, sd, sesd), qnorm(0.95, sd, sesd) ) )

})
