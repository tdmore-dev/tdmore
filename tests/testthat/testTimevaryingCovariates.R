library(tdmore)
library(testthat)
library(RxODE)
library(tdmore)
library(magrittr)
library(dplyr)

context("Test the time-varying covariates")

modelCode <- "
ALB_Obs = ALB;
K=0.05;
d/dt(centr) = -K*centr;
"
rxModel <- RxODE::RxODE(modelCode)
ev <- eventTable()
ev$add.dosing(dose=5*72, start.time = 0)
ev$add.sampling(time=seq(0, 42))
rxModel$solve(params=c(ALB=42), events=ev)

### WARNING: Time-varying covariates should use the time defined by the eventTable's ev$get.EventTable
## Execute using RxODE with time-varying covariates: covariates frame not according to spec (not enough time points)
covariates <- data.frame(time=c(0, 2, 20), ALB=c(100, 80, 40))
ev <- eventTable()
ev$add.dosing(dose=5*72, start.time = 0)
ev$add.sampling(time=seq(0, 42))
rxModel$solve(params=c(ETAK=0, ETAV=0), events=ev, covs = covariates, covs_interpolation="constant")
## Conclusion: RxODE gives the wrong prediction!

## Execute using RxODE, sampling times equal to number of observations
covariates <- data.frame(time=seq(0, 42, by=1), ALB=100*runif(43))
ev <- eventTable()
ev$add.dosing(dose=5*72, start.time = 0)
ev$add.sampling(time=seq(0, 42))
out <- rxModel$solve(params=c(ETAK=0, ETAV=0), events=ev, covs = covariates, covs_interpolation="constant") %>% as.data.frame
out %>% right_join(covariates, by="time")
## Right prediction

## Execute using RxODE, sampling times equal to number of observations
covariates <- data.frame(time=seq(0, 42, by=0.5), ALB=100*runif(85))
ev <- eventTable()
ev$add.dosing(dose=5*72, start.time = 0)
ev$add.sampling(time=seq(0, 42))
out <- rxModel$solve(params=c(ETAK=0, ETAV=0), events=ev, covs = covariates, covs_interpolation="constant") %>% as.data.frame
out %>% right_join(covariates, by="time") %>% mutate(ALB_calc = ALB/42)
## Right prediction

## Conclusion: Using covariates requires some rewriting of the input data.frame

modelCode <- "
TVK_1 = 0.0512;
TVK_2 = 0.0540;
TVK_3 = 0.0667;
COV_CRP_K = 0.0831;
COV_SA_K = -0.825;
TVV_BASE = 6.87;
COV_SEX_V = 1.30;

if(SEX == 0) {
  TVV = TVV_BASE;
} else {
  TVV = TVV_BASE*COV_SEX_V;
}
if(MPRE == 1) {
  TVK_BASE = TVK_1;
} else if (MPRE == 2) {
  TVK_BASE = TVK_2;
} else {
  TVK_BASE = TVK_3;
}

TVK = TVK_BASE*(ALB/42.0)**(COV_SA_K)*(CRP/6.1)**(COV_CRP_K);

K = TVK * exp(ETAK*sqrt(0.0856));
V = TVV * exp(ETAV*sqrt(0.0852));

CONC = centr / V;
d/dt(centr) = -K*centr;
d/dt(AUC) = CONC;
"
rxModel <- RxODE::RxODE(modelCode)
res_var <- list(errorModel(var = "CONC", add=sqrt(0.09), prop=sqrt(0.0376)))
model <- rxModel %>% tdmore(parameters=c("ETAK", "ETAV"), res_var, covs_interpolation="constant")
regimen <- data.frame(
  TIME=c(0, 14, 42),
  AMT=5*72 #5mg
)
observed <- data.frame(TIME=10, CONC=20)
covariates <- data.frame(TIME=c(0, 5), ALB=c(30, 42), CRP=c(4, 6), SEX=0, MPRE=2)

#debugonce(tdmore:::model_predict.RxODE)
#debugonce(tdmore:::predict.tdmore)
#debugonce(zoo::na.locf)
model %>% predict(newdata=data.frame(TIME=84, CONC=NA), regimen=regimen, params=c(ETAK=0, ETAV=0), covariates=covariates)


pred <- model %>% estimate(regimen=regimen, covariates=covariates)
ipred <- model %>% estimate(observed=observed, regimen=regimen, covariates=covariates)

library(ggplot2)
ggplot(observed, aes(x=TIME, y=CONC)) + geom_point() +
  geom_line(aes(color="Population"), data=predict(pred, newdata=seq(0, 16, by=0.1))) +
  geom_line(aes(color="Individual"), data=predict(ipred, newdata=seq(0, 16, by=0.1)))
