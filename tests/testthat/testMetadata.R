library(tdmore)
library(nlmixr)
library(testthat)
library(magrittr)


context("Test the TDMore metadata")

output1 <- output(name="CONC", label="Drug concentration 1", unit="ng/ml", default_value=5)
output2 <- output(name="CONC2", label="Drug concentration 2", unit="ng/ml")

cov1 <- covariate(name="BW", label="Weight", unit="kg", min=30, max=80)
cov2 <- covariate(name="CYP3A5", label="CYP3A5 expressor", choices=list(Fast=0, Slow=1))

dose <- dose(unit="mg", dosing_interval=12, default_value=5)

target <- target(min=10, max=15)

test_that("To string methods work as expected", {
  expect_equal(toString(output1), "Drug concentration 1 (ng/ml)")
  expect_equal(toString(cov1), "Weight (kg)")
  expect_equal(toString(dose), "mg")
  expect_equal(toString(target), "Target: [10,15]")
})

test_that("Metadata are compatible with TDMore", {
  m1 <- (meropenem_model_wt) %>% tdmore() %>% metadata(output1, output2) %>% metadata(cov1, cov2) %>% metadata(dose) %>% metadata(target)

  expect_equal(m1$metadata, list(output1, output2, cov1, cov2, dose, target))
  expect_equal(getMetadataByName(m1, "CONC"), output1)
  expect_equal(getMetadataByName(m1, "DOSE"), dose)
  expect_equal(getMetadataByName(m1, "TARGET"), target)
})


