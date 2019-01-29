library(tdmore)
library(nlmixr)
library(ggplot2)
library(testthat)

context("Test that the mixture models are working")

set.seed(0)

regimen <- data.frame(
  TIME=c(0,24),
  AMT=150
)

mod_1cpt_1 <- nlmixrUI(function(){
  ini({
    TVV <- 70
    TVKA <- 1
    TVCL <- 4
    ECL ~ 0.09 # SD=0.3
    EPS_PROP <- 0.1
  })
  model({
    CL <- TVCL * exp(ECL)
    V <- TVV
    KA <- TVKA
    d/dt(abs) = - abs*KA
    d/dt(center) = abs*KA - CL/V * center
    CONC = center / V
    CONC ~ prop(EPS_PROP)
  })
})
m1 <- (mod_1cpt_1) %>% tdmore()

mod_1cpt_2 <- nlmixrUI(function(){
  ini({
    TVV <- 70
    TVKA <- 1
    TVCL <- 2
    ECL ~ 0.09 # SD=0.3
    EPS_PROP <- 0.1
  })
  model({
    CL <- TVCL * exp(ECL)
    V <- TVV
    KA <- TVKA
    d/dt(abs) = - abs*KA
    d/dt(center) = abs*KA - CL/V * center
    CONC = center / V
    CONC ~ prop(EPS_PROP)
  })
})
m2 <- (mod_1cpt_2) %>% tdmore()

plot(m1, regimen)
plot(m2, regimen)

expect_error(tdmore_mixture(m1, m2, probs = c(0.5, 0.51))) # Error: Sum of probabilities must be 1
mixture <- tdmore_mixture(m1, m2, probs = c(0.5, 0.5))

# Model 1 is the most likely model
observed1 <- data.frame(TIME=c(10, 20), CONC=c(1.20, 0.75))
ipred_mixture <-  estimate(object=mixture, observed = observed1, regimen = regimen)
plot(ipred_mixture)
expect_true(ipred_mixture$winner == 1)
expect_equal(round(ipred_mixture$fits_prob$IPk, digits = 4), c(0.9248, 0.0752))

# Model 2 is the most likely model
observed2 <- data.frame(TIME=c(10, 20), CONC=c(1.5, 1.25))
ipred_mixture <-  estimate(object=mixture, observed = observed2, regimen = regimen)
plot(ipred_mixture)
expect_true(ipred_mixture$winner == 2)
expect_equal(round(ipred_mixture$fits_prob$IPk, digits = 4), c(0.1500, 0.8500))

# Model 2 is the most likely model
observed3 <- data.frame(TIME=c(10, 15, 20), CONC=c(1.5, 1.5, 0.75))
ipred_mixture <-  estimate(object=mixture, observed = observed3, regimen = regimen)
plot(ipred_mixture)
expect_true(ipred_mixture$winner == 2)
expect_equal(round(ipred_mixture$fits_prob$IPk, digits = 4), c(0.4910, 0.5090))

# Model 2 was the most likely but probabilities are now different, so model 1 is the most likely model
mixture <- tdmore_mixture(m1, m2, probs = c(0.6, 0.4))
observed3 <- data.frame(TIME=c(10, 15, 20), CONC=c(1.5, 1.5, 0.75))
ipred_mixture <-  estimate(object=mixture, observed = observed3, regimen = regimen)
plot(ipred_mixture)
expect_true(ipred_mixture$winner == 1)
expect_equal(round(ipred_mixture$fits_prob$IPk, digits = 4), c(0.5913, 0.4087))

# It should be also possible to predict from a mixture, default model index is used
data1 <- predict(
  object = mixture,
  newdata = data.frame(TIME = seq(0, 24, by = 0.5), CONC = NA),
  regimen = regimen,
  se = TRUE
)

#---------------------------------------------------
# Test integration of mixture models with tdmore set
#---------------------------------------------------

m3 <- (meropenem_model_wt) %>% tdmore()
set <- tdmore_set(m3, mixture)

# Plot should be model 1
ipred <-  estimate(object=set, observed = observed1, regimen = regimen)
plot(ipred)

# Plot should be model 3 because covariates in m3 model
ipred <-  estimate(object=set, observed = observed1, regimen = regimen, covariates = c(WT=70))
plot(ipred)
