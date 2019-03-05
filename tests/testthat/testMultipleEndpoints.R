library(tdmore)
library(nlmixr)
library(testthat)

context("Test the Sunitinib PK/PD model")
set.seed(0)

# Load the default tdmore
tdmore <- (sunitinib_pkpd_model) %>% tdmore(maxsteps=1E3*500)

# Checking the error model
expect_output(print(tdmore$res_var[[1]]) ,
              "CONC  : proportional  \\( prop= 0\\.417  \\)")
expect_output(print(tdmore$res_var[[2]]),
              "SLD.*proportional.*prop= 0\\.143.*"
              )

# Estimation example
regimen <- data.frame(
  TIME=0,
  AMT=50,
  II=24,
  ADDL=40*7
)
times <- seq(0, 40*7*24, by=1)
predict(tdmore, regimen=regimen, newdata=seq(0, 30*7*24, by=7*24))
observed <- data.frame(TIME=c(0, 30*7*24), CONC=NA, SLD=c(25, 14))
ipred <- tdmore %>% estimate(observed = observed, regimen = regimen)
expect_equal(round(coef(ipred)['EBASE'], digits=4) , c(EBASE=0.5397))
