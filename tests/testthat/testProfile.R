library(tdmore)
library(RxODE)
library(nlmixr)
library(ggplot2)
library(testthat)
library(dplyr)

context("Test that the profile method works as intended")

profile <- function(...) {
  set.seed(0)
  stats::profile(...)
}

# Load the fluticasone model with 5 IIV parameters
tdmore <- getModel("fluticasone")

# Predicting new data
regimen <- data.frame(
  TIME=c(0, 8, 16),            # Every 8 hour and for 1 day, an injection is given
  AMT=c(1000, 1000, 1000),     # 1g is administered
  RATE=c(1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
)

# Estimating individual parameters
pred <- tdmore %>% estimate(regimen = regimen)
observed <- data.frame(TIME=c(9, 16), CONC=c(55, 40))

ipred <- tdmore %>% estimate(observed = observed, regimen = regimen)
plot(ipred)

# Compute profile on 3 parameters
profile <- profile(ipred, maxpts = 10, fix=coef(ipred)[c("EV2", "EQ")])

expect_error( plot(profile), regexp="Cannot produce plot in more than 2 dimensions")
expect_warning({
  z <- plot(profile, parameters=c("ECL", "EV1"))
}, regexp="multiple .*scores for a given parameter combination")

z <- z + facet_wrap(~EKa)
#print(z) #Ka does not shift the LL-surface of ECL/EV1
vdiffr::expect_doppelganger("profile_with_Ka_facet", z)

z <- plot(profile, parameters="EKa")
#print(z) #Ka clearly does not matter
vdiffr::expect_doppelganger("profile_on_Ka", z)

# You can specify NA to estimate the parameter at every step
profile <- profile(ipred, maxpts=20, fix=c(EKa=0, coef(ipred)[c("EV2", "EQ")]))
z1 <- plot(profile, parameters="ECL") #LL-profile with EV1 jumping around
vdiffr::expect_doppelganger("ecl_different_ev1", z1)

profile <- profile(ipred, maxpts=20, fix=c(EV1=NA, EKa=0, coef(ipred)[c("EV2", "EQ")]))
z2 <- plot(profile) #LL-profile with EV1 estimated
vdiffr::expect_doppelganger("ecl_optimized_ev1", z2)

#gridExtra::grid.arrange(z1 + coord_cartesian(xlim=c(-0.5, 0.5), ylim=c(0, 0.003)),
#                        z2 + coord_cartesian(xlim=c(-0.5, 0.5), ylim=c(0, 0.003)))

expect_error({
  plot(profile, parameters=list() )
})

# Error if requesting a parameter that does not exist
expect_error(
  plot(profile, parameters="NO_EXIST")
)

profile <- profile(ipred, maxpts = 20, limits = list(ECL=c(-0.3,0.3), EV1=c(-0.8,0)), fix=coef(ipred)[c("EKa", "EV2", "EQ")])
vdiffr::expect_doppelganger("profile_limits", plot(profile, contour=F))
