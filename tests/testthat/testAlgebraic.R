
myFunction <- function(t, A0, A1, TIME, AMT, EKA, EV, ECL, WT) {
  V = 70 * exp(EV) * WT/70
  CL = 4 * exp(ECL) * (WT/70)^0.75
  KA = 1 * exp(EKA)
  K = CL / V
  tD = TIME

  pk1cptoral()
}
m1 <- algebraic(myFunction, inits=c(A0=0, A1=0))
regimen <- data.frame(TIME=seq(0, 100, by=24), AMT=150)
covariates <- c(WT=49)

test_that("model_predict function generates values as expected", {
  expect_error( #parameter missing
    model_predict(model=m1, times=numeric(), parameters=c(EKA=0, EV=0), covariates=covariates)
  )
  #one parameter too many: is ignored
  model_predict(model=m1, times=numeric(), parameters=c(EKA=0, EV=0, ECL=0, EV3=0), covariates=covariates, regimen=regimen)

  expect_error( #typo in a parameter
    model_predict(model=m1, times=numeric(), parameters=c(EKA=0, EV1=0, ECL=0), covariates=covariates)
  )
  expect_error( #covariates missing
    model_predict(model=m1, times=numeric(), parameters=c(EKA=0, EV=0, ECL=0))
  )
  expect_equal(
    model_predict(model=m1, times=numeric(), parameters=c(EKA=0, EV=0, ECL=0),
                  regimen=regimen, covariates=covariates),
    data.frame(TIME=numeric(), A0=numeric(), A1=numeric(), CONC=numeric())
  )
  prediction <- model_predict(model=m1, times=0:15,
                              parameters=c(EKA=0, EV=0, ECL=0), covariates=covariates,
                              regimen=regimen)
  expect_equal(
    prediction %>% colnames,
    c("TIME", "A0", "A1", "CONC")
  )
  expect_equal(prediction$TIME, 0:15)
  expect_equal(prediction$CONC, myFunction(t=0:15, A0=0, A1=0, TIME=0, AMT=150, EKA=0, EV=0, ECL=0, WT=49)$CONC )

  #prediction <- model_predict(model=m1, times=0:15,
  #                parameters=c(EKA=0, EV=0, ECL=0), covariates=covariates,
  #                regimen=data.frame(TIME=c(0, 12), AMT=c(NA, 1000)))
  #plot(prediction, type='l')

  result <- model_predict(model=m1, times=seq(0, 100, by=0.1),
                          regimen=regimen,
                          parameters=c(EKA=0, EV=0, ECL=0),
                          covariates=covariates)
  z1 <- ggplot(result, aes(x=TIME, y=CONC)) + geom_line()
  expect_doppelganger("Algebraic function model_predict", z1)
})

m2 <- tdmore(m1,
             res_var=errorModel(prop=0.1),
             omega=c(EKA=0.3, EV=0.3, ECL=0.3))
test_that("predict.tdmore() function generates values as expected", {
  result <- predict(m2, newdata=data.frame(TIME=seq(0, 100, by=0.1), CONC=NA), regimen=regimen, covariates=covariates)
  z1 <- ggplot(result, aes(x=TIME, y=CONC)) + geom_line()
  expect_doppelganger("Algebraic function predict", z1)
})




myFunction <- function(t, TIME, AMT, DURATION, EV, ECL, WT) {
  V = 70 * exp(EV) * WT/70
  CL = 4 * exp(ECL) * (WT/70)^0.75
  K = CL / V
  RATE = AMT / DURATION
  pk1cptinfusion()
}
m1 <- algebraic(myFunction)
regimen <- data.frame(TIME=seq(0, 100, by=24), AMT=150)
covariates <- c(WT=49)

test_that("More exotic errors with extra treatment aspects", {
  expect_error({
    model_predict(model=m1, times=0:15,
                              parameters=c(EV=0, ECL=0), covariates=covariates,
                              regimen=regimen)
  }, regexp=".DURATION. is missing")

  prediction <- model_predict(model=m1, times=0:15,
                parameters=c(EV=0, ECL=0), covariates=covariates,
                regimen=data.frame(TIME=999, AMT=100, DURATION=300))
  expect_equal(prediction$CONC, rep(0, 16))

  prediction <- model_predict(model=m1, times=0:15,
                              parameters=c(EV=0, ECL=0), covariates=covariates,
                              regimen=data.frame(TIME=3, AMT=100, DURATION=300))
  expect_equal(prediction$CONC, c(0, 0, 0, 0, 0.00659458697316676, 0.0127897997728064, 0.0186098248614116,
                                  0.0240773839473502, 0.0292138226917053, 0.0340391940429475, 0.0385723365247844,
                                  0.0428309477828265, 0.0468316536772008, 0.0505900731908523, 0.0541208794069401,
                                  0.0574378567933874))
})

myFunction <- function(t, TIME, AMT, DURATION, EKA, EV, ECL, WT) {
  V = 70 * exp(EV) * WT/70
  CL = 4 * exp(ECL) * (WT/70)^0.75
  KA = 1 * exp(EKA)
  K = CL / V
  tD = TIME
  pk1cptoral()
}
m1 <- algebraic(myFunction)
test_that("IOV works", {
  prediction <- model_predict(model=m1, times=0:15,
                              parameters=c(EKA=-10, EKA=10, EV=0, ECL=0), covariates=covariates,
                              regimen=data.frame(TIME=c(0, 8), AMT=100, DURATION=300, OCC=c(1,2)), iov=c("EKA"))
  expect_error({
    model_predict(model=m1, times=0:15,
              parameters=c(EKA=-10, EKA=10, EV=0, ECL=0), covariates=covariates,
              regimen=data.frame(TIME=c(0, 8, 12), AMT=100, DURATION=300, OCC=c(1,2, 3)), iov=c("EKA"))
  }) #not enough parameters for 3 occasions

  model_predict(model=m1, times=0:15,
              parameters=c(EKA=0, EV=0, ECL=0), covariates=data.frame(TIME=c(0, 8), WT=c(300, 40)),
              regimen=data.frame(TIME=c(0, 8), AMT=100, DURATION=300))
  expect_snapshot_value(prediction, style="serialize")
})
