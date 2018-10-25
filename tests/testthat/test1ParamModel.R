library(tdmore)
library(RxODE)
library(nlmixr)
library(ggplot2)
library(testthat)

context("Test a 1-parameter-model")

set.seed(0)

# Load the Meropenem model
source(paste0(test_path(), ("/modelLibrary.R")))
tdmore <- nlmixrUI(meropenem_1param_model) %>% tdmore()

# Predicting new data
regimen <- data.frame(
  TIME=c(0, 8, 16),            # Every 8 hour and for 1 day, an injection is given
  AMT=c(1000, 1000, 1000),     # 1g is administered
  RATE=c(1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
)

# Estimate and plot IPRED
observed <- data.frame(TIME=c(9, 16), CONC=c(30, 8))
ipred <- estimate(tdmore = tdmore, observed = observed, regimen = regimen)

plot(ipred)
plot(profile(ipred))

expect_equal(round(ipred$res, digits=3), c(EV1=0.578))
