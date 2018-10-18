library(tdmore)
library(RxODE)
library(nlmixr)
library(ggplot2)
library(testthat)
library(dplyr)

context("Test that the profile method works as intended")

set.seed(0)

# Creating your model

modelCode <- function(){
  ini({
    TVV1 <- 24.4;
    TVV2 <- 7.01;
    TVQ <- 4.97;
    TVCL <- 9.87;
    ECL ~ 0.194 # This value corresponds to OMEGA_CL (44% SD)
    EV1 ~ 0.287 # This value corresponds to OMEGA_V1 (54% SD)
    EPS_PROP <- 0.371 # Proportional error (37% SD)
  })
  model({
    CL <- TVCL * exp(ECL)
    V1 <- TVV1 * exp(EV1)
    V2 <- TVV2
    Q <- TVQ
    K12 <- Q/V1
    K21 <- Q/V2

    d/dt(center) = - CL/V1 * center - K12*center + K21 * periph
    d/dt(periph) = K12*center - K21 * periph

    CONC = center / V1
    CONC ~ prop(EPS_PROP) # Proportional error linked to the PK model
  })
}

nlmixrUI <- nlmixrUI(modelCode)
tdmore <- tdmore(nlmixrUI)

# Predicting new data

regimen <- data.frame(
  TIME=c(0, 8, 16),            # Every 8 hour and for 1 day, an injection is given
  AMT=c(1000, 1000, 1000),     # 1g is administered
  RATE=c(1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
)

# Estimating individual parameters

pred <- estimate(tdmore = tdmore, regimen = regimen)
observed <- data.frame(TIME=c(9, 16), CONC=c(30, 6))

ipred <- estimate(tdmore = tdmore, observed = observed, regimen = regimen)
plot(ipred)

# Compute different profiles

profile <- profile(ipred, maxpts = 20)
plot(profile)
plot(profile, raster = F)
plot(profile, contour = F)

profile <- profile(ipred, maxpts = 20, limits = list(ECL=c(-0.8,0.5), EV1=c(-0.8,0.9)))
plot(profile)

profile <- profile(ipred, maxpts = 20, limits = list(ECL=c(-0.5,0)))
plot(profile)

profile <- profile(ipred, maxpts = 20, limits = list(EV1=c(-1,1)), fix=c(ECL=-0.15))
plot(profile)
