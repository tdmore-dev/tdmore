library(tdmore)
library(RxODE)
library(nlmixr)
library(ggplot2)
library(testthat)

context("Test that the nlmixr class works as intended")

modelCode <- function(){
  ini({
    TVKA <- 3.7
    TVV1 <- 61
    TVQ <- 10
    TVCL <- 3.7
    ECL ~ 0.0784 #ETA1, 28%
    EV1 ~ 0.0361 #ETA2, 19%
    EPS_PROP <- 0.23 # Proportional error, 23% SD
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

# Create the observed and covariates dataframe
observed <- data.frame(TIME=c(2), CONC=c(0.040))
covariates <- data.frame(TIME=c(0), WT=c(70))

# Compute PRED
pred <- model %>% estimate(regimen = regimen, covariates = covariates)
expect_equal(pred$res, c(ECL=0.0, EV1=0.0))

# Compute PRED
pred <- model %>% estimate(regimen = regimen, covariates = covariates)
expect_equal(pred$res, c(ECL=0.0, EV1=0.0))

# Compute IPRED
ipred <- model %>% estimate(observed = observed, regimen = regimen, covariates = covariates)
expect_equal(round(ipred$res, digits=4), c(ECL=0.0336, EV1=0.1175))

# Custom plot, compare PRED & IPRED
ggplot(observed, aes(x=TIME, y=CONC)) + geom_point() +
  geom_line(aes(color="Population"), data=predict(pred, newdata=seq(0, 48, by=0.1))) +
  geom_line(aes(color="Individual"), data=predict(ipred, newdata=seq(0, 48, by=0.1)))

# Plot IPRED
#debugonce(tdmore:::plot.tdmorefit)
p1 <- plot(ipred, newdata=data.frame(TIME=seq(0, 48, by=0.1), CONC=NA))


# Find dose example and test
observed <- data.frame(TIME=c(2,26), CONC=c(0.040, 0.0675))
covariates <- data.frame(TIME=c(0), WT=c(70))

# Compute new IPRED
ipred <- model %>% estimate(observed = observed, regimen = regimen, covariates = covariates)

# Plot IPRED
p2 <- plot(ipred, newdata=data.frame(TIME=seq(0, 48, by=0.1), CONC=NA))

# Find dose test
D <- findDose(ipred, regimen=regimen, target=data.frame(TIME=35, CONC=0.0395))
expect_equal(round(D$root, digits=2), 5.84)

plot(p1)
plot(p2)
