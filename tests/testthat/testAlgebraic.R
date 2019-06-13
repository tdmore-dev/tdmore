context("Test algebraic models")

myFunction <- function(t, TIME, AMT, EKA, EV, ECL, WT) {
  V = 70 * exp(EV) * WT/70
  CL = 4 * exp(ECL) * (WT/70)^0.75
  KA = 1 * exp(EKA)
  k = CL / V
  tD = TIME
  ifelse(t >= TIME,
         AMT/V * (KA/(KA-k)) * (exp(-k*(t-tD)) - exp(-KA*(t-tD))),
         0)
}
m1 <- algebraic(myFunction)
regimen <- data.frame(TIME=seq(0, 100, by=24), AMT=150)
covariates <- c(WT=49)

test_that("model_predict function generates values as expected", {
  expect_error( #parameter missing
    model_predict(model=m1, times=numeric(), parameters=c(EKA=0, EV=0), covariates=covariates)
  )
  expect_error( #one parameter too many
    model_predict(model=m1, times=numeric(), parameters=c(EKA=0, EV=0, ECL=0, EV3=0), covariates=covariates, regimen=regimen)
  )
  expect_error( #typo in a parameter
    model_predict(model=m1, times=numeric(), parameters=c(EKA=0, EV1=0, ECL=0), covariates=covariates)
  )
  expect_error( #covariates missing
    model_predict(model=m1, times=numeric(), parameters=c(EKA=0, EV=0, ECL=0))
  )
  expect_equal(
    model_predict(model=m1, times=numeric(), parameters=c(EKA=0, EV=0, ECL=0),
                  regimen=regimen, covariates=covariates),
    data.frame(TIME=numeric(), CONC=numeric())
  )
  prediction <- model_predict(model=m1, times=0:15,
                              parameters=c(EKA=0, EV=0, ECL=0), covariates=covariates,
                              regimen=regimen)
  expect_equal(
    prediction %>% colnames,
    c("TIME", "CONC")
  )
  expect_equal(prediction$TIME, 0:15)
  expect_equal(prediction$CONC, myFunction(t=0:15, TIME=0, AMT=150, EKA=0, EV=0, ECL=0, WT=49) )

  expect_error(
    model_predict(model=m1, times=0:15,
                  parameters=c(EKA=0, EV=0, ECL=0), covariates=covariates,
                  regimen=data.frame(TIME=c(0, 12), AMT=c(NA, 1000))),
    "The provided regimen contains NA.*"
  )

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




myFunction <- function(t, TIME, AMT, DURATION, EKA, EV, ECL, WT) {
  V = 70 * exp(EV) * WT/70
  CL = 4 * exp(ECL) * (WT/70)^0.75
  KA = 1 * exp(EKA)
  k = CL / V
  tD = TIME
  ifelse(t >= TIME,
         AMT/V * (KA/(KA-k)) * (exp(-k*(t-tD)) - exp(-KA*(t-tD))),
         0)
}
m1 <- algebraic(myFunction)
regimen <- data.frame(TIME=seq(0, 100, by=24), AMT=150)
covariates <- c(WT=49)

test_that("More exotic errors with extra treatment aspects", {
  expect_error({
    model_predict(model=m1, times=0:15,
                              parameters=c(EKA=0, EV=0, ECL=0), covariates=covariates,
                              regimen=regimen)
  })

  prediction <- model_predict(model=m1, times=0:15,
                parameters=c(EKA=0, EV=0, ECL=0), covariates=covariates,
                regimen=data.frame(TIME=Inf, AMT=100, DURATION=300))
  expect_equal(prediction$CONC, rep(0, 16))
})

myFunction <- function(t, TIME, AMT, DURATION, EKA, EV, ECL, WT) {
  V = 70 * exp(EV) * WT/70
  CL = 4 * exp(ECL) * (WT/70)^0.75
  KA = 1 * exp(EKA)
  k = CL / V
  tD = TIME
  ifelse(t >= TIME,
         AMT/V * (KA/(KA-k)) * (exp(-k*(t-tD)) - exp(-KA*(t-tD))),
         0)
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
  })

  prediction <- model_predict(model=m1, times=0:15,
                parameters=c(EKA=0, EV=0, ECL=0), covariates=data.frame(TIME=c(0, 8), WT=c(300, 40)),
                regimen=data.frame(TIME=c(0, 8), AMT=100, DURATION=300))
  expect_known_value(prediction, "algebraic-different-weights.txt")
})


myFunction <- function(t, TIME, AMT, EKA, EV, ECL, WT) {
  V = 70 * exp(EV) * WT/70
  CL = 4 * exp(ECL) * (WT/70)^0.75
  KA = 1 * exp(EKA)
  k = CL / V
  tD = TIME
  CONC <- ifelse(t >= TIME,
         AMT/V * (KA/(KA-k)) * (exp(-k*(t-tD)) - exp(-KA*(t-tD))),
         0)
  list(CONC=CONC, V=V, CL=CL, KA=KA)
}
m1 <- algebraic(myFunction, output="CONC")
test_that("multiple model outputs", {
  prediction <- model_predict(model=m1, times=0:15,
                              parameters=c(EKA=-10, EKA=10, EV=0, ECL=0), covariates=covariates,
                              regimen=data.frame(TIME=c(0, 8), AMT=100, OCC=c(1,2)), iov="EKA")
  expect_known_value(prediction, "algebraic_with_multiple_outputs.txt")
})
