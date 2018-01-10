## Goal: test defining structural model using RxODE
library(testthat)
library(tdmore)
context("Test that the Model class works as intended")

library(RxODE)
library(plyr)
library(dplyr)
library(ggplot2)

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
d/dt(centr) = ka*abs - k12*perip + k21*perip - ke*centr;
d/dt(perip) = k12*perip - k21*perip;
"
model <- RxODE::RxODE(modelCode) %>%
  tdmore(prop=0.23) #Model has 23% proportional error

regimen <- data.frame(
  TIME=seq(0, 7)*24,
  AMT=5 #5mg
)
observed <- data.frame(TIME=2, CONC=0.04)

pred <-  model %>%
  estimate(regimen=regimen)
ipred <- model %>%
  estimate(observed, regimen)

print(summary(ipred))

profile <- profile(fit, maxpts=20)
ggplot(profile, aes(x=ETA1, y=ETA2, z=logLik)) + geom_contour()

newdata = data.frame(TIME=seq(0, 12, length.out=50), CONC=NA)

ggplot(fit %>% predict(newdata), aes(x=TIME, y=CONC)) +
  geom_line(aes(color="Fit")) +
  geom_ribbon(aes(fill="Fit (95% CI)", ymin=CONC.lower, ymax=CONC.upper), data=fit %>% predict(newdata, se.fit=TRUE), alpha=0.3)+
  geom_line(aes(color="Population"), data=pred %>% predict(newdata)) +
  geom_point(aes(color="Observed"), data=observed)




