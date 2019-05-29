library(tidyverse)
library(tdmore)
library(nlmixr)
library(testthat)

context("Test that MPC estimation works as intended")

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
covariates <- as.data.frame(t(c(TIME=0, theta, WT=70)))

regimen <- data.frame(
  TIME=c(0, 24, 48, 72, 96),
  AMT=5,
  OCC=c(1,2,2,3,3)
)

observed <- data.frame(
  TIME=c(24, 48, 72, 96, 120),
  CONC=c(2.5, 4, 5, 5, 4.5)
)

plot(m1, newdata=seq(0, 5*24, by=0.1), regimen=regimen, covariates=covariates, se.fit=NA)

# Estimate with EBE
# plot(estimate(m1, regimen=regimen, covariates=covariates, observed=observed))

# Estimate with MPC
m1_mpc <- m1 %>% mpc(theta=theta, suffix="_next")

#debugonce(tdmore:::estimate.tdmore_mpc)
#debug(tdmore:::model_prepare.RxODE)
ipred <- estimate(m1_mpc, regimen=regimen, covariates=covariates, observed=observed)

expectedCoef <- c(ECL=-0.2018, EV1=-0.0559, ECL=-0.1534, EV1=-0.0331, ECL=0.1174, EV1=-0.0016)
expect_equal(round(coef(ipred), digits=4), expectedCoef)

# PLOT
times <- seq(0, max(observed$TIME), by=0.2)
fit <- predict(m1, newdata=times, regimen=regimen, parameters=coef(ipred), covariates=ipred$covariates)
pred <- predict(m1, newdata=times, regimen=regimen, parameters=c(), covariates=covariates)

plot <- ggplot(mapping=aes(x=TIME, y=CONC)) +
  geom_line(data=pred, color="blue") +
  geom_line(data=fit, color="red") +
  geom_point(data=observed)

print(plot)

# Need a special function to plot ipred mcp?
class(ipred$tdmore) <- "tdmore" # TEMPORARY
plot(ipred, se.fit=F)


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
fit <- predict(m1, newdata=times, regimen=regimen, parameters=coef(ipred), covariates=ipred$covariates)
pred <- predict(m1, newdata=times, regimen=regimen, parameters=c(), covariates=covariates)

plot <- ggplot(mapping=aes(x=TIME, y=CONC)) +
  geom_line(data=pred, color="blue") +
  geom_line(data=fit, color="red") +
  geom_point(data=observed)

print(plot)
