library(tdmore)
library(nlmixr)
library(testthat)
library(magrittr)


output1 <- output(name="CONC", label="Drug concentration 1", unit="ng/ml", default_value=5)
output2 <- output(name="CONC2", label="Drug concentration 2", unit="ng/ml")

cov1 <- covariate(name="BW", label="Weight", unit="kg", min=30, max=80)
cov2 <- covariate(name="CYP3A5", label="CYP3A5 expressor", choices=c(Fast=0, Slow=1))

formulation1 <- formulation(name="Compound",unit="mg", dosing_interval=12, default_value=5)
formulation2 <- formulation(name="Compound2",unit="mg", dosing_interval=8, default_value=6)

target <- target(min=10, max=15)

test_that("To string methods work as expected", {
  expect_output(print(output1), "Drug concentration 1 \\(ng/ml\\)")
  expect_output(print(cov1), "Weight \\(kg\\)")
  expect_output(print(formulation1), "mg")
  expect_output(print(target), "Target: \\[10,15\\]")
})

test_that("Metadata are compatible with TDMore", {
  m1 <- getModel("example_wt") %>% metadata(output1, output2) %>% metadata(cov1, cov2) %>% metadata(formulation1,formulation2) %>% metadata(target)

  expect_equal(m1$metadata, list(CONC=output1, CONC2=output2, BW=cov1, CYP3A5=cov2, Compound=formulation1, Compound2=formulation2, TARGET=target))
  expect_equal(getMetadataByName(m1, "CONC"), output1)
  expect_equal(getMetadataByName(m1, getMetadataByClass(m1,"tdmore_formulation", all=T)[[1]]$name), formulation1)
  expect_equal(getMetadataByName(m1, getMetadataByClass(m1,"tdmore_formulation", all=T)[[2]]$name), formulation2)
  expect_equal(getMetadataByName(m1, "TARGET"), target)
})


