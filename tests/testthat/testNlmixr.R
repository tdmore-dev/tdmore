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
    OV1 <- 0.28
    OCL <- 0.19
    ECL ~ 1 # ETA1
    EV1 ~ 1 #ETA2
    EPS_PROP <- 0.23
  })
  model({
    KA <- TVKA
    CL <- TVCL * (70/70)^0.75 * exp(ECL*OCL)
    V1 <- TVV1 * exp(EV1*OV1)
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

model <-nlmixrUI(modelCode) %>% tdmore()

regimen <- data.frame(
  TIME=seq(0, 7)*24,
  AMT=5 #5mg
)
observed <- data.frame(TIME=2, CONC=0.04)

pred <- model %>% estimate(regimen=regimen)
ipred <- model %>% estimate(observed, regimen)

pred %>% predict()
ipred %>% predict()

newdata = data.frame(TIME=seq(0, 12, length.out=50), CONC=NA)

