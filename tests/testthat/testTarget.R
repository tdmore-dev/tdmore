#getTroughs <- function(model, regimen, deltamin=1/4, deltaplus=1/4)
myModel <- algebraic(fun=function(t, TIME, AMT, EV) {pk1cptiv_(t=t, TIME=TIME, AMT=AMT, V=10, K=0.2)}) %>%
  tdmore(omega=c(EV=0.1), res_var=errorModel(add=0.1)) %>%
  metadata(
    formulation("CompoundA", dosing_interval=12, unit="mg")
  )

describe("getTroughs shifts time to ensure trough, not peak", {
  myModel2 <- RxODE::RxODE("
  K=0.2;
  V=10 * exp(EV);
  d/dt(A1) = -K*A1;
  CONC=A1/V;
  ") %>% tdmore::tdmore(omega=c(EV=0.1), res_var=list(tdmore::errorModel("CONC", prop=0.1))) %>%
    tdmore::metadata( formulation("CompoundA", dosing_interval=12, unit="mg") )

  regimen <- tibble::tibble(TIME=c(0, 12), AMT=10, FORM="CompoundA")
  profile <- predict(myModel2, regimen=regimen, newdata=seq(0, 13, by=0.01))
  trough <- predict(myModel2, regimen=regimen, newdata=12-1e-15)
  expect_lt(trough$CONC, min(profile$CONC)) #trough is lowest concentration!

  peak <- predict(myModel2, regimen=regimen, newdata=12)
  expect_gt(peak$CONC, max(profile$CONC)) #observation at t=12 for IV profile is MAX concentration

  trough <- predict(myModel2, regimen=regimen, newdata=getTroughs(myModel2, regimen)[1] )
  expect_lt(trough$CONC, min(profile$CONC)) #trough is lowest concentration!
})

describe("getTroughs", {
  it("requires a FORM column", {
    expect_error(getTroughs(myModel, tibble::tibble(
      TIME=c(24, 48)
      )))
  })
  it("requires all FORM to exist as metadata", {
    expect_error(getTroughs(myModel, tibble::tibble(
      TIME=c(24, 48),
      FORM="Unknown"
    )), "Formulation.*not defined.*")
  })

  it("calculates the right supposed troughs in the future", {
    expect_equal(getTroughs(myModel, tibble::tibble(
      TIME=c(24),
      FORM="CompoundA"
    )), 24+12)
  })

  it("does not match another treatment if too far off", {
    expect_equal(getTroughs(myModel, tibble::tibble(
      TIME=c(24, 100),
      FORM="CompoundA"
    )), c(24, 100)+12)
  })

  it("matches another treatment if in bounds", {
    expect_equal(getTroughs(myModel, tibble::tibble(
      TIME=c(24, 24+10),
      FORM="CompoundA"
    )), c(24+10, 24+10+12))
  })

  it("warning if another treatment too early", {
    expect_warning(x <- getTroughs(myModel, tibble::tibble(
      TIME=c(24, 24+6),
      FORM="CompoundA"
    )), "A treatment was planned.*")
    expect_equal(x, c(30, 42))
  })
})
