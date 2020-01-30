#getTroughs <- function(model, regimen, deltamin=1/4, deltaplus=1/4)
myModel <- algebraic(fun=function(t, TIME, AMT, EV) {pk1cptiv_(t=t, TIME=TIME, AMT=AMT, V=10, K=0.2)}) %>%
  tdmore(omega=c(EV=0.1), res_var=errorModel(add=0.1)) %>%
  metadata(
    formulation("CompoundA", dosing_interval=12, unit="mg")
  )

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

  it("errors if another treatment too early", {
    expect_error(getTroughs(myModel, tibble::tibble(
      TIME=c(24, 24+6),
      FORM="CompoundA"
    )), "A treatment was detected .*")
  })
})
