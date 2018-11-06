library(tdmore)
library(nlmixr)
library(testthat)

context("Test the Sunitinib PK/PD model")
set.seed(0)

# Load the default tdmore
source(paste0(test_path(), ("/modelLibrary.R")))
tdmore <- nlmixrUI(sunitinib_pkpd_model) %>% tdmore(maxsteps=1E3*500)

# Checking the error model
expect_equal(unlist(tdmore$res_var[1]) , c(var="CONC", add=0, prop=0.417, exp=0))
expect_equal(unlist(tdmore$res_var[2]) , c(var="SLD", add=0, prop=0.143, exp=0))

# Estimation example
regimen <- data.frame(
  TIME=0,
  AMT=50,
  II=24,
  ADDL=40*7
)
times <- seq(0, 40*7*24, by=1)
observed <- data.frame(TIME=c(0, 30*7*24), CONC=NA, SLD=c(25, 14))
ipred <- estimate(tdmore = tdmore, observed = observed, regimen = regimen)
expect_equal(round(coef(ipred)['EBASE'], digits=4) , c(EBASE=0.5397))

