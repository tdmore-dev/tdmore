library(tdmore)
library(nlmixr)
library(ggplot2)
library(testthat)

# Load the default tdmore
tdmore <- getModel("default")

regimen <- data.frame(
  TIME=seq(0, 1)*24,
  AMT=5 #5mg
)

expect_error(plot(tdmore, regimen, vars=c("CONC"))) # Error: Value for covariate(s) WT missing
expect_error(plot(tdmore, regimen, vars=c("CONC"), covariates=c(WT=50)), NA) # No error is raised
expect_error(tdmore %>% predict(newdata=0:24, regimen=regimen)) # Error: Value for covariate(s) WT missing
expect_error(tdmore %>% predict(newdata=0:24, regimen=regimen, covariates=c(WT=50)), NA) # No error is raised
