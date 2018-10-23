library(tdmore)
library(RxODE)
library(nlmixr)
library(ggplot2)
library(testthat)

context("When OMEGA is 0, then the ETA should always be fixed to 0 as well.")

# Load the Meropenem model
source(paste0(test_path(), ("/modelLibrary.R")))
tdmore <- nlmixrUI(meropenem_omega0_model) %>% tdmore(covs_interpolation="constant")

# Predicting new data
regimen <- data.frame(
  TIME=c(0, 8, 16),            # Every 8 hour and for 1 day, an injection is given
  AMT=c(1000, 1000, 1000),     # 1g is administered
  RATE=c(1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
)

observed <- data.frame(TIME=c(9, 16), CONC=c(30, 8))

ipred <- estimate(tdmore = tdmore, observed = observed, regimen = regimen)

plot(ipred)

profile <- profile(ipred, maxpts = 50)
plot(profile)
