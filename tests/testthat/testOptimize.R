m1 <- nlmixrUI(function(){
  ini({
    TVCL <- 3.7
    TVV1 <- 61
    ECL ~ 0.0784 #ETA1, 28%
    EV1 ~ 0.0361 #ETA2, 19%
    EPS_PROP <- 0.23 # Proportional error, 23% SD
  })
  model({
    KA <- 3.7
    CL <- TVCL * exp(ECL)
    V1 <- TVV1 * exp(EV1)
    V2 <- V1
    Q <- 10
    K12 <- Q/V1
    K21 <- Q/V2

    d/dt(depot) = -KA*depot
    d/dt(center) = KA*depot - CL/V1 * center - K12*center + K21 * periph
    d/dt(periph) = K12*center - K21 * periph

    CONC = (center / V1) * 100
    CONC ~ prop(EPS_PROP)
  })
}) %>% tdmore(iov=c("EV1", "ECL")) %>%
  metadata(formulation("CompoundA", dosing_interval=24, unit="mg"),
           target(min=5, max=10))

regimen <- data.frame(
  TIME=c(0, 24, 48, 72, 96, 120),
  AMT=5,
  OCC=c(1,2,2,3,3,4),
  FORM="CompoundA"
)
observed <- data.frame(
  TIME=c(24, 48, 72, 96, 120)-0.5,
  CONC=c(2.5, 4, 5, 5, 4.5)
)

ipred <- estimate(m1, regimen=regimen, observed=observed)
plot(ipred)
describe("findDoses", {
  it("finds the right dose to hit a target", {
    myRegimen <- regimen
    myRegimen$FIX <- TRUE
    myRegimen$FIX[2] <- FALSE
    rec <- findDoses(ipred, regimen=myRegimen)

    newRegimen <- myRegimen
    newRegimen$AMT[2] <- 14.009638
    expect_equal(rec$regimen, newRegimen, tolerance=1e-2)
  })
  it("is equivalent to findDose for a single target", {
    myRegimen <- regimen
    myRegimen$FIX <- TRUE
    myRegimen$FIX[2] <- FALSE
    rec1 <- findDoses(ipred, regimen=myRegimen)
    target <- tibble::tibble(TIME=modifyMantissa(myRegimen$TIME[2]+24, a=-2^-52), CONC=7.5)

    rec2 <- findDose(ipred, regimen=myRegimen, doseRows=which(!myRegimen$FIX), target=target)
    expect_equal(rec1, rec2)
  })
  it("Optimizes the whole treatment", {
    rec <- findDoses(ipred)
    res <- predict(ipred, regimen=rec$regimen)
    expect_equal(res$CONC, rep(7.5, 5), tolerance=0.1) #equal to 7.5, more or less

    plot(ipred) + autolayer(rec)
  })
})
