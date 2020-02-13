library(tdmore)
library(nlmixr)
library(testthat)
library(dplyr)
library(ggplot2)

context("Test that MPC estimation works as intended")

describe("addIterationColumn", {
  it("adds the right iteration number", {
    regimen <- data.frame(TIME=seq(0, 7*24, by=24), AMT=50, OCC=1:8)
    observed <- data.frame(TIME=c(23, 45, 47, 6*24-1) )
    expect_equal(
      tdmore:::addIterationColumn(regimen, observed)$ITER,
      c(1,2,2,3)
  )})

  it("pad 1 to the left part of the observations", {
    regimen <- data.frame(TIME=50+seq(0, 7*24, by=24), AMT=50, OCC=1:8)
    observed <- data.frame(TIME=c(23, 45, 47, 6*24-1) )
    expect_equal(
      tdmore:::addIterationColumn(regimen, observed)$ITER,
      c(1,1,1,2)
    )
  })
})


m1 <- nlmixrUI(function(){
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
}) %>% tdmore(iov=c("EV1", "ECL"))

theta <- c(TVCL=3.7, TVV1=61)

covariatesThetaIncluded <- as.data.frame(t(c(TIME=0, theta, WT=70)))

regimen <- data.frame(
  TIME=c(0, 24, 48, 72, 96, 120),
  AMT=5,
  OCC=c(1,2,2,3,3,4)
)

observed <- data.frame(
  TIME=c(24, 48, 72, 96, 120)-0.5,
  CONC=c(2.5, 4, 5, 5, 4.5)
)


expect_equal(
  tdmore:::addIterationColumn(regimen, observed),
  cbind(observed, ITER=c(1,2,2,3,3))
)

plot(m1, newdata=seq(0, 5*24, by=1), regimen=regimen, covariates=covariatesThetaIncluded) +
  geom_point(data=observed, aes(x=TIME, y=CONC))

# Estimate with EBE
ebe <- estimate(m1, regimen=regimen, covariates=covariatesThetaIncluded, observed=observed)
plot(ebe, newdata=0:150) + geom_label(data=regimen, aes(x=TIME, y=0, label=OCC), hjust=-1)
coef(ebe)
# last coefficients should be 0
expect_equal(
  tail(coef(ebe), n=2),
  c(ECL=0, EV1=0),
  tol=1e-6
)

# Estimate with MPC
m1_mpc <- m1 %>% mpc(theta=theta, suffix="_next")

#debugonce(tdmore:::estimate.tdmore_mpc)

covariates <- data.frame(TIME=0, WT=70) # No need to give the MPC thetas for an MPC model
ipred <- estimate(m1_mpc, regimen=regimen, covariates=covariates, observed=observed)

expectedCoef <- c(ECL = -0.191879467130269, EV1 = -0.0545878097352673,
                  ECL = -0.15158796212222,  EV1 = -0.0330037469432886,
                  ECL = 0.116116940276883, EV1 = -0.00132757305153023,
                  ECL = 0, EV1 = 0) #as many coef as occasions!
expect_equal(coef(ipred), expectedCoef, tol=1e-4)

# PLOT
times <- seq(0, max(observed$TIME), by=0.2)
fit <- predict(m1_mpc, newdata=times, regimen=regimen, parameters=coef(ipred), covariates=ipred$covariates)
pred <- predict(m1_mpc, newdata=times, regimen=regimen, parameters=c(), covariates=covariates)

plot <- ggplot(mapping=aes(x=TIME, y=CONC)) +
  geom_line(data=pred, color="blue") +
  geom_line(data=fit, color="red") +
  geom_point(data=observed)

print(plot)

plot(ipred, se.fit=FALSE, newdata=0:150)
plot(ipred, se.fit=TRUE, newdata=0:150)

# ---------------------------------------------------------------------

regimen <- data.frame(
  TIME=c(0, 24),
  AMT=5,
  OCC=c(1,2)
)

observed <- data.frame(
  TIME=c(13),
  CONC=c(5)
)

ipred <- estimate(m1_mpc, regimen=regimen, covariates=covariates, observed=observed)

times <- seq(0, 48, by=0.2)
fit <- predict(m1_mpc, newdata=times, regimen=regimen, parameters=coef(ipred), covariates=ipred$covariates)
pred <- predict(m1_mpc, newdata=times, regimen=regimen, parameters=c(), covariates=covariates)

plot <- ggplot(mapping=aes(x=TIME, y=CONC)) +
  geom_line(data=pred, color="blue") +
  geom_line(data=fit, color="red") +
  geom_point(data=observed)

print(plot)


# ---------------------------------------------------------------------

# Time-varying covariates

timeVaryingCovs <- data.frame(TIME=c(0,20,40), WT=c(50,70,90))

# No observation
ipred1 <- estimate(m1_mpc, regimen=regimen, covariates=timeVaryingCovs)
expect_equal(ipred1$covariates, data.frame(TVCL=c(3.7, 3.7, 3.7),
                                           TVV1=c(61, 61, 61),
                                           TIME=c(0, 20, 40),
                                           WT=c(50,70,90)))

# Observation
ipred2 <- estimate(m1_mpc, regimen=regimen, covariates=timeVaryingCovs, observed=observed)
expect_true(all.equal(ipred2$covariates, data.frame(TVCL=c(3.7, 3.7, 2.7878, 2.7878),
                                           TVV1=c(61, 61,  47.8952,  47.8952),
                                           TIME=c(0, 20, 24, 40),
                                           WT=c(50, 70, 70, 90)), check.attributes=F, tolerance=1E-4))

times <- seq(0, 48, by=0.2)
fit <- predict(m1_mpc, newdata=times, regimen=regimen, parameters=coef(ipred2), covariates=ipred2$covariates)
pred <- predict(m1_mpc, newdata=times, regimen=regimen, parameters=c(), covariates=timeVaryingCovs)

plot <- ggplot(mapping=aes(x=TIME, y=CONC)) +
  geom_line(data=pred, color="blue") +
  geom_line(data=fit, color="red") +
  geom_point(data=observed)

print(plot)



