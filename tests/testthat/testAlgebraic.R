
context("Test your own function")

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
result <- model_predict(model=m1, newdata=data.frame(TIME=seq(0, 100, by=0.1), CONC=NA),
              regimen=regimen,
              parameters=c(EKA=0, EV=0, ECL=0))
z1 <- ggplot(result, aes(x=TIME, y=CONC)) + geom_line()
vdiffr::expect_doppelganger("Algebraic function manual model_predict plot", z1)

m2 <- tdmore(m1,
             res_var=errorModel(prop=0.1),
             omega=c(EKA=0.3, EV=0.3, ECL=0.3))
predict(m2, newdata=data.frame(TIME=seq(0, 100, by=0.1), CONC=NA),
        regimen=regimen)
z1 <- plot(m2, newdata=data.frame(TIME=seq(0, 100), CONC=NA), regimen=regimen)
vdiffr::expect_doppelganger("Algebraic function tdmore plot", z1)
