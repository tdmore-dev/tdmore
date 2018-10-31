library(tdmore)
library(RxODE)
library(nlmixr)
library(ggplot2)
library(testthat)
library(dplyr)

context("Test that the profile method works as intended")

set.seed(0)

# Load the Meropenem model
source(paste0(test_path(), ("/modelLibrary.R")))
tdmore <- nlmixrUI(meropenem_model) %>% tdmore()

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
