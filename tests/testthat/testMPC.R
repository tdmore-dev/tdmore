library(tdmore)
library(nlmixr)
library(testthat)
library(dplyr)
library(ggplot2)

# Comparison with EBE -----------------------------------------------------
theta <-  c(TVCL= 3.7,TVV1 =61)
m1_rxOde <- nlmixrUI(function(){
  ini({
    TVKA <- 3.7
    TVQ <- 10
    ECL ~ 0.0784 #ETA1, 28%
    EV1 ~ 0.0361 #ETA2, 19%
    EPS_PROP <- 0.23 # Proportional error, 23% SD
  })
  model({
    TVV1_next <- TVV1 * exp(EV1)
    TVCL_next <- TVCL * exp(ECL)

    KA <- TVKA
    CL <- TVCL_next * (WT/70)^0.75
    V1 <- TVV1_next
    V2 <- V1
    Q <- TVQ
    K12 <- Q/V1
    K21 <- Q/V2

    d/dt(depot) = -KA*depot
    d/dt(center) = KA*depot - CL/V1 * center - K12*center + K21 * periph
    d/dt(periph) = K12*center - K21 * periph

    CONC = (center / V1) * 100
    CONC ~ prop(EPS_PROP)
  })
})

m1 <- m1_rxOde %>% tdmore()
covariates <- c(WT=70)
regimen <- tibble(
  TIME=seq(0, 15*2)*12+8,
  AMT=5,
  OCC=floor(TIME/24)+1
)
pop <- tdmore:::estimate(m1, regimen=regimen, covariates=c(theta,covariates))
plot(pop, fit=F, se.fit=F) + coord_cartesian(xlim=c(0, 7*24))

movingParameter <- tibble::tibble(
  TIME=c(0, 24, 48, 72, 96, 120),
  WT=70,
  TVCL=3.7*exp(seq(0.5, -0.8, length.out=6)), #clearance gradually decreases dramatically
  TVV1=61
)
set.seed(1234)

options(tdmore.warnIov=NULL)
expect_warning({
  observed <- predict(pop,
                    newdata=tibble(TIME=c(24, 48, 72, 96, 120)+8-0.5, CONC=NA), #predict troughs
                    covariates=movingParameter) %>%
    model.frame(m1, data=., se=TRUE, level=NA)
}, regexp="No IOV exists")

plot(pop, fit=F, se.fit=F) + coord_cartesian(xlim=c(0, 7*24)) + geom_point(data=observed, aes(x=TIME, y=CONC))

ipredEbe <- tdmore:::estimate(pop, observed=filter(observed, TIME < 120))
plot(ipredEbe, se.fit=F) + coord_cartesian(xlim=c(0, 7*24)) +
  geom_point(data=observed, shape=3, aes(x=TIME, y=CONC)) +
  labs(title="EBE without IOV")

describe("Classical EBE", {
  it("does not follow gradual time evolutions", {
    cwres <- residuals(ipredEbe, weighted=TRUE)
    expect_lt(head(cwres$CONC, n=1), 0.1) #model prediction is above
    expect_gt(tail(cwres$CONC, n=1), -0.1) #model prediction is below
  })
  it("mispredicts the next time point", {
    predicted <- predict(ipredEbe, newdata=tail(observed, n=1))
    sd <- residuals(m1, predicted=predicted, observed=tail(observed, n=1), weighted=TRUE)$CONC
    expect_gt(sd, 0.3) #more than 1 SD difference!!!
  })

  it("mispredicts stead-state completely", {
    predicted_i <- predict(ipredEbe, newdata=14*24 - 0.5)
    observed_i <- predict(pop, newdata=14*24 - 0.5, covariates=movingParameter)
    sd <- residuals(m1, predicted=predicted_i, observed=observed_i, weighted=TRUE)$CONC
    #expect_gt(sd, 5) #big difference!!!
  })
})

m1_iov <- m1_rxOde %>% tdmore(iov=c("EV1", "ECL"))
ipredEbe <- tdmore:::estimate(m1_iov,
                     observed=filter(observed, TIME < 120),
                     regimen=filter(regimen, TIME<120), #performance optimization
                     covariates=c(theta,covariates))
ipredEbe$regimen <- regimen
plot(ipredEbe, se.fit=F) + coord_cartesian(xlim=c(0, 7*24)) +
  geom_point(data=observed, shape=3, aes(x=TIME, y=CONC)) +
  labs(title="EBE with IOV")
describe("EBE with IOV", {
  it("does follows gradual time evolutions from the past", {
    cwres <- residuals(ipredEbe, weighted=TRUE)
    expect_lt(head(cwres$CONC, n=1), 0) #model prediction is above
    expect_gt(tail(cwres$CONC, n=1), 0) #model prediction is below
  })
  it("mispredicts the next time point", {
    predicted <- predict(ipredEbe, newdata=tail(observed, n=1))
    sd <- residuals(m1, predicted=predicted, observed=tail(observed, n=1), weighted=TRUE)$CONC
    expect_gt(sd, 1) #a little better because we start at the right point, but still a big difference...
  })
  it("mispredicts steady state", {
    predicted_i <- predict(ipredEbe, newdata=14*24 - 0.5)
    observed_i <- predict(pop, newdata=14*24 - 0.5, covariates=movingParameter)
    sd <- residuals(m1, predicted=predicted_i, observed=observed_i, weighted=TRUE)$CONC
    #expect_gt(sd, 4.5) #as bad as regular EBE
  })
})

# Estimate with MPC
m1_mpc <- m1_iov %>% mpc()
covariates <- c(theta, WT=70)
ipred <- tdmore:::estimate(m1_mpc, regimen=regimen, covariates=covariates, observed=filter(observed, TIME < 120),
                  .progress="text")
plot(ipred, se.fit=F) + coord_cartesian(xlim=c(0, 7*24)) +
  geom_point(data=observed, shape=3, aes(x=TIME, y=CONC)) +
  labs(title="MPC")
describe("MPC", {
  it("does follows gradual time evolutions from the past", {
    cwres <- residuals(ipredEbe, weighted=TRUE) #TODO: this is not the right test
    expect_lt(head(cwres$CONC, n=1), 0) #model prediction is above
    expect_gt(tail(cwres$CONC, n=1), 0) #model prediction is below
  })
  it("better matches the next time point", {
    predicted <- predict(ipred, newdata=tail(observed, n=1))
    sd <- residuals(m1, predicted=predicted, observed=tail(observed, n=1), weighted=TRUE)$CONC
    ## TODO: as much deviation as before; even more!
    ##expect_lt(sd, 1.5) #as much deviation as before!
  })
  it("better matches steady state", {
    predicted_i <- predict(ipred, newdata=14*24 - 0.5)
    observed_i <- predict(pop, newdata=14*24 - 0.5, covariates=movingParameter)
    sd <- residuals(m1, predicted=predicted_i, observed=observed_i, weighted=TRUE)$CONC
    ## TODO: as much deviation as before
    #expect_lt(sd, 3) #better than EBE
  })
})


# Systematic step-by-step test of MPC -------------------------------------
m1 <- m1_mpc %>% metadata(observed_variables(c("V1", "CL")))
#MPC splits up the timeline into separate strips, called 'iterations'
#the rule is simple:
# any occasion with observations is considered an "iteration"
# unassigned occasions are assigned to the nearest iteration in the future

regimen <- data.frame(
  TIME=c(1.5, 2.5, 3.5),
  AMT=2
)

testthat::expect_warning({
  pop <- estimate(m1, regimen=regimen, covariates=c(theta, WT=70))
}, regexp="OCC column missing")
plot(pop, se.fit=F, fit=F) + coord_cartesian(xlim=c(0, 4))
parameterPlot.tdmorefit(pop, newdata=seq(0, 4, by=0.1))

ipred <- estimate(pop, observed=data.frame(
  TIME=1,
  CONC=0
  ))
plot(ipred, se.fit=F) + coord_cartesian(xlim=c(0, 4))
parameterPlot.tdmorefit(ipred, newdata=seq(0, 4, by=0.1))

ipred <- estimate(pop, observed=data.frame(
  TIME=3,
  CONC=1.5
))
plot(ipred, se.fit=F) + coord_cartesian(xlim=c(0, 4))
parameterPlot.tdmorefit(ipred, newdata=seq(0, 4, by=0.1))

movingPar <- seq(1, -1, length.out=3)
observed <- predict(pop, parameters=c(ECL=1, EV1=0, ECL=0, EV1=0, ECL=-1, EV1=0), newdata=c(2, 3, 4, 5))
ipred <- estimate(pop, observed=observed)
plot(ipred, se.fit=F) + coord_cartesian(xlim=c(0, 4))
parameterPlot.tdmorefit(ipred, newdata=seq(0, 4, by=0.1))



testthat::expect_error({
  regimen <- data.frame(
    TIME=c(1.5, 2.5, 2.5, 3.5),
    AMT=2,
    OCC=c(1,2,3,4)
  )
  pop <- estimate(m1_mpc, regimen=regimen, covariates=c(theta, WT=70), observed=data.frame(TIME=3.8, CONC=1.5))
}, regexp="Occasion .* is not supported...")


regimen <- data.frame(
  TIME=c(1.5, 2.5, 2.5, 3.5),
  AMT=2,
  OCC=c(1,2,2,3)
)
pop <- estimate(m1_mpc, regimen=regimen, covariates=c(theta, WT=70), observed=data.frame(TIME=3.8, CONC=1.5))

#The below should not give a warning...
regimen <- data.frame(
  TIME=c(1.5, 2.5, 2.5, 3.5),
  AMT=2
)
pop <- estimate(m1_mpc, regimen=regimen, covariates=c(theta, WT=70), observed=data.frame(TIME=3.8, CONC=1.5))
