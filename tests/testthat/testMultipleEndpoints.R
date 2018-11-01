library(tdmore)
library(nlmixr)
library(ggplot2)
library(testthat)

context("Test multiple endpoints models")

# Load the default tdmore
source(paste0(test_path(), ("/modelLibrary.R")))
#debugonce(tdmore:::tdmore.nlmixrUI)
tdmore <- nlmixrUI(two_error_models) %>% tdmore()

regimen <- data.frame(
  TIME=seq(0, 1)*24,
  AMT=5 #5mg
)
covariates = c(WT=70)

# Default tdmore plot
set.seed(1)
plot(tdmore, regimen, covariates=covariates, mc.maxpts = 100)
