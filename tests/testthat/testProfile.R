library(tdmore)
library(RxODE)
library(nlmixr)
library(ggplot2)
library(testthat)
library(dplyr)

context("Test that the profile method works as intended")

set.seed(0)

# Load the Meropenem model
tdmore <- (meropenem_model) %>% tdmore()

# Predicting new data
regimen <- data.frame(
  TIME=c(0, 8, 16),            # Every 8 hour and for 1 day, an injection is given
  AMT=c(1000, 1000, 1000),     # 1g is administered
  RATE=c(1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
)

# Estimating individual parameters
pred <- tdmore %>% estimate(regimen = regimen)
observed <- data.frame(TIME=c(9, 16), CONC=c(30, 6))

ipred <- tdmore %>% estimate(observed = observed, regimen = regimen)
plot(ipred)

# Compute different profiles
# TODO: Use vdiffr to check values
profile <- profile(ipred, maxpts = 20)
plot(profile)
plot(profile, raster = F)
plot(profile, contour = F)
plot(profile, parameters="ECL") #TODO: Gives the wrong figure

# You can specify NA to estimate the parameter at every step
profile <- profile(ipred, maxpts = 20, fix=c(ECL=NA))
plot(profile)
# Or simply provide a FIX parameter so it is the same throughout
profile <- profile(ipred, maxpts = 20, fix=coef(ipred)['ECL'])
plot(profile)

#TODO: specifying 3 parameters does not give an error message
#Instead, the first 2 are drawn while silently ignoring the 3rd parameter
plot(profile, parameters=c("ECL", "EV1") )
expect_warning({
  plot(profile, parameters=list() )
})
# Error if requesting a parameter that does not exist
expect_error(
  plot(profile, parameters="NO_EXIST")
)

profile <- profile(ipred, maxpts = 20, limits = list(ECL=c(-0.8,0.5), EV1=c(-0.8,0.9)))
plot(profile)

profile <- profile(ipred, maxpts = 20, limits = list(ECL=c(-0.5,0)))
plot(profile)

profile <- profile(ipred, maxpts = 20, limits = list(EV1=c(-1,1)), fix=c(ECL=-0.15))
plot(profile)
