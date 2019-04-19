library(tdmore)
library(nlmixr)
library(testthat)
library(magrittr)


context("Test the TDMore metadata")

output1 <- output("CONC", label="Drug concentration 1", unit="ng/ml")
output2 <- output("CONC2", label="Drug concentration 2", unit="ng/ml")

cov1 <- covariate("BW", label="Weight", unit="kg", min=30, max=80)
cov2 <- covariate("CYP3A5", label="CYP3A5 expressor", choices=list(Fast=0, Slow=1))

m1 <- (meropenem_model_wt) %>% tdmore() %>% metadata(output1, output2) %>% metadata(cov1, cov2)

expect_equal(m1$metadata, list(output1, output2, cov1, cov2))
expect_equal(getMetadataByName(m1, "CONC"), output1)
