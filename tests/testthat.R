Sys.setenv("R_TESTS" = "")
library(testthat)
library(vdiffr)
library(tdmore)

test_check("tdmore")
