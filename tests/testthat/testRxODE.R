## Goal: test defining structural model using RxODE
rm(list=ls(all=T))
library(testthat)
context("Test that the Model class works as intended")

library(RxODE)
## TODO: this does not compile correctly in R CMD CHECK, NOOOOOO :-(
modelCode <- "
CL = 3.7 * exp(ETA1*0.19);
Vc = 61 * exp(ETA2*0.28);
ka=3.7;
Q = 10;
Vp = Vc;
k12=Q/Vc;
k21=Q/Vp;
ke=CL/Vc;

CONC = centr / Vc;

d/dt(abs) = -ka*abs;
d/dt(centr) = ka*abs - k12*centr + k21*perip - ke*centr;
d/dt(perip) = k12*centr - k21*perip;
"

model <- RxODE::RxODE(modelCode) %>%
  tdmore(prop=0.23) #Model has 23% proportional error

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
