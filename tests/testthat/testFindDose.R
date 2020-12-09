library(tdmore)
library(RxODE)
library(nlmixr)
library(ggplot2)
library(testthat)

set.seed(0)

describe("findDose works correctly", {
  # Load the example model
  tdmore <- getModel("example")
  # Predicting new data
  regimen <- data.frame(
    TIME=c(0, 8, 16),            # Every 8 hour and for 1 day, an injection is given
    AMT=c(1000, 1000, 1000),     # 1g is administered
    RATE=c(1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
  )
  data <- predict(
    object = tdmore,
    newdata = data.frame(TIME = seq(0, 24, by = 0.5), CONC = NA),
    regimen = regimen,
    se = TRUE
  )
  pred <- tdmore %>% estimate(regimen = regimen)
  observed <- data.frame(TIME=c(9, 16), CONC=c(30, 8))
  ipred <- tdmore %>% estimate(observed = observed, regimen = regimen)
  data <- predict(ipred, newdata=data.frame(TIME=seq(0, 24, 0.1), CONC=NA), se=TRUE)

  #autoplot(ipred) + coord_cartesian(xlim=c(0, 25))
  # Finding the right dose to give
  newRegimen <- data.frame(
    TIME=c(0, 8, 16, 24),              # A fourth dose on the second day is added
    AMT=c(1000, 1000, 1000, NA),       # Adding an unknown dose on the second day
    RATE=c(1000, 1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
  )

  it("does not allow findDose on multiple targets", {
    testthat::expect_error({
      findDose(
        ipred,
        regimen = newRegimen,
        interval = c(100, 5000),
        target = data.frame(TIME = c(32, 32+24), CONC = c(3.1, 4.1))
      )
    }, regexp="multiple targets")

    testthat::expect_error({
      findDose(
        ipred,
        regimen = newRegimen,
        interval = c(100, 5000),
        target = data.frame(TIME = 32, CONC = 3.1, periph=20)
      )
    }, regexp="multiple targets")
  })

  it("gives correct dose recommendations", {
    recommendation1 <- findDose(
      ipred,
      regimen = newRegimen,
      interval = c(100, 5000),
      target = data.frame(TIME = 32, CONC = 3.1)
    )
    expect_equal(round(as.double(recommendation1), digits=2), 296.62)
  })

  it("allows getting an SE on dose", {
    # Finding the right dose to give (with 95% CI)
    recommendation2 <- suppressWarnings(findDose(
      ipred,
      regimen = newRegimen,
      interval = c(100, 5000),
      target = data.frame(TIME = 32, CONC = 3.1),
      se.fit = T,
      level = 0.95,
      mc.maxpts = 50
    ))

    expect_equal(unlist(round(recommendation2$dose, digits=1)), c(
      dose.median = 224.8,
      dose.lower = 100,
      dose.upper = 775.8
    ), tolerance=0.5) ## high tolerance, even new RxODE versions result in different values...
  })

  it("updateRegimen works correctly", {
    modRegimen <- updateRegimen( regimen, newDose=1234)
    expect_equal(modRegimen$AMT, c(1000, 1000, 1234) )
  })
})
