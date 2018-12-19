
context("Test algebraic models")

m1 <- algebraic({
  function(t, TIME, AMT, EKA, EV, ECL) {
    V = 70 * exp(EV)
    CL = 4 * exp(ECL)
    KA = 1 * exp(EKA)
    k = CL / V
    tD = TIME
    ifelse(t >= TIME,
           AMT/V * (KA/(KA-k)) * (exp(-k*(t-tD)) - exp(-KA*(t-tD))),
           0)
  }
})
regimen <- data.frame(TIME=seq(0, 100, by=24), AMT=150)
test_that("model_predict function generates values as expected", {
  expect_equal(
    model_predict(model=m1, times=numeric(), parameters=c(EKA=0, EV=0, ECL=0)),
    data.frame(TIME=numeric(), CONC=numeric())
  )
  expect_equal(
    model_predict(model=m1, times=0:15, parameters=c(EKA=0, EV=0, ECL=0)),
    data.frame(TIME=0:15, CONC=0)
  )

  result <- model_predict(model=m1, times=seq(0, 100, by=0.1),
                          regimen=regimen,
                          parameters=c(EKA=0, EV=0, ECL=0))
  z1 <- ggplot(result, aes(x=TIME, y=CONC)) + geom_line()
  vdiffr::expect_doppelganger("Algebraic function model_predict", z1)
})

m2 <- tdmore(m1,
             res_var=errorModel(prop=0.1),
             omega=c(EKA=0.3, EV=0.3, ECL=0.3))
test_that("predict.tdmore() function generates values as expected", {
  result <- predict(m2, newdata=data.frame(TIME=seq(0, 100, by=0.1), CONC=NA), regimen=regimen)
  z1 <- ggplot(result, aes(x=TIME, y=CONC)) + geom_line()
  vdiffr::expect_doppelganger("Algebraic function predict", z1)
})

