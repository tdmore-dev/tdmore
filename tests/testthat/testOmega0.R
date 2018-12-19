library(tdmore)
library(RxODE)
library(nlmixr)
library(ggplot2)
library(testthat)

context("When OMEGA is 0, an error is raised")

# Load Meropenem model with omega 0
expect_error((meropenem_omega0_model) %>% tdmore()) # Error is raised

# Load default meroponem model (no omega 0)
expect_error((meropenem_1param_model) %>% tdmore(), NA) # No error is raised
