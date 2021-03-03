library(tdmore)
library(nlmixr)
library(testthat)

set.seed(0)

# Load the default tdmore
m1 <- getModel("example_pkpd") %>% tdmore(maxsteps=1E3*500)

# Checking the error model
expect_output(print(m1$res_var[[1]]) ,
              "CONC  : proportional  \\( prop= 0\\.4  \\)")
expect_output(print(m1$res_var[[2]]),
              "SLD.*proportional.*prop= 0\\.15"
              )

# Estimation example
regimen <- data.frame(
  TIME=0,
  AMT=50,
  II=24,
  ADDL=40*7
)
times <- seq(0, 40*7*24, by=1)
predict(as.population(m1), regimen=regimen, newdata=seq(0, 30*7*24, by=7*24))
observed <- data.frame(TIME=c(0, 30*7*24), CONC=NA, SLD=c(25, 14))
ipred <- m1 %>% tdmore:::estimate(observed = observed, regimen = regimen)
expect_equal(round(coef(ipred)['EBASE'], digits=4) , c(EBASE=0.5687))
