Sys.setenv("R_TESTS" = "")
library(testthat)
library(tdmore)

test_check("tdmore")
