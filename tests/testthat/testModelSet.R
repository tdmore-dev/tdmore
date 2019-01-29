library(tdmore)
library(RxODE)
library(nlmixr)
library(ggplot2)
library(testthat)

context("Test that the model set class is working")

set.seed(0)

# Load the Meropenem model with WT covariate
m1 <- (meropenem_model_wt) %>% tdmore()

# Load the Meropenem model without covariates
m2 <- (meropenem_model) %>% tdmore()

# Create a model set
set <- tdmore_set(m1, m2)

# Predicting new data
regimen <- data.frame(
  TIME=c(0, 8, 16),            # Every 8 hour and for 1 day, an injection is given
  AMT=c(1000, 1000, 1000),     # 1g is administered
  RATE=c(1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
)

# Predict without covariates
data1 <- predict(
  object = set,
  newdata = data.frame(TIME = seq(0, 24, by = 0.5), CONC = NA),
  regimen = regimen,
  se = TRUE
)

ggplot(data1, aes(x=TIME, y=CONC)) +
  geom_ribbon(aes(fill="Population", ymin=CONC.lower, ymax=CONC.upper), fill="steelblue2", alpha=0.15) +
  geom_line(aes(color="Population"), data=data1) +
  scale_color_manual(values=c("steelblue2"))

pred <- estimate(object = set, regimen = regimen)

# Check tdmore model in pred does not have covariates
expect_true(length(pred$tdmore$covariates)==0)

# Predict with covariates
data2 <- predict(
  object = set,
  newdata = data.frame(TIME = seq(0, 24, by = 0.5), CONC = NA),
  regimen = regimen,
  covariates = c(WT=25),
  se = TRUE
)

ggplot(data2, aes(x=TIME, y=CONC)) +
  geom_ribbon(aes(fill="Population", ymin=CONC.lower, ymax=CONC.upper), fill="steelblue2", alpha=0.15) +
  geom_line(aes(color="Population"), data=data2) +
  scale_color_manual(values=c("steelblue2"))

pred <- estimate(object = set, regimen = regimen, covariates = c(WT=25))

# Check tdmore model in pred does have covariates
expect_true(pred$tdmore$covariates==c("WT"))


