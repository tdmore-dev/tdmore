## Goal: test defining structural model using nlmixr
rm(list=ls(all=T))
library(testthat)
context("Test that the Model class works as intended")

library(nlmixr)

modelCode <- function(){
  ini({
    TVKA <- 3.7
    TVV1 <- 61
    TVQ <- 10
    TVCL <- 3.7
    ECL ~ 0.0784 #ETA1 (0.28^2)
    EV1 ~ 0.0361 #ETA2 (0.19^2)
    EPS_PROP <- 0.23
  })
  model({
    KA <- TVKA
    CL <- TVCL * (WT/70)^0.75 * exp(ECL)
    V1 <- TVV1 * exp(EV1)
    V2 <- V1
    Q <- TVQ
    K12 <- Q/V1
    K21 <- Q/V2

    d/dt(depot) = -KA*depot
    d/dt(center) = KA*depot - CL/V1 * center - K12*center + K21 * periph
    d/dt(periph) = K12*center - K21 * periph

    CONC = center / V1
    CONC ~ prop(EPS_PROP)
  })
}

model <- nlmixrUI(modelCode) %>% tdmore(covs_interpolation="constant")

regimen <- data.frame(
  TIME=seq(0, 1)*24,
  AMT=5 #5mg
)

observed <- data.frame(TIME=c(2), CONC=c(0.040))
covariates <- data.frame(TIME=c(0), WT=c(70))

# Compute PRED
pred <- model %>% estimate(regimen = regimen, covariates = covariates)
stopifnot(all.equal(pred$res, c(ECL=0.0, EV1=0.0)))

# Compute IPRED
ipred <- model %>% estimate(observed = observed, regimen = regimen, covariates = covariates)
stopifnot(all.equal(round(ipred$res, digits=4), c(ECL=0.0336, EV1=0.1175)))

# Custom plot, compare PRED & IPRED
library(ggplot2)
ggplot(observed, aes(x=TIME, y=CONC)) + geom_point() +
  geom_line(aes(color="Population"), data=predict(pred, newdata=seq(0, 48, by=0.1))) +
  geom_line(aes(color="Individual"), data=predict(ipred, newdata=seq(0, 48, by=0.1)))

# Plot IPRED
plot(ipred, newdata=data.frame(TIME=seq(0, 48, by=0.1), CONC=NA))

#----------------------------------------------------------------------

observed <- data.frame(TIME=c(2,26), CONC=c(0.040, 0.0675))
covariates <- data.frame(TIME=c(0), WT=c(70))

# Compute new IPRED
ipred <- model %>% estimate(observed = observed, regimen = regimen, covariates = covariates)

# Plot IPRED
plot(ipred, newdata=data.frame(TIME=seq(0, 48, by=0.1), CONC=NA))

# Find dose test
D <- findDose(ipred, regimen=regimen, target=data.frame(TIME=35, CONC=0.0395))
