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
model <- RxODE::RxODE(modelCode)


regimen <- data.frame(
  TIME=seq(0, 7)*24,
  AMT=5 #5mg
)
observed <- data.frame(TIME=2, CONC=0.04)

pred <- model %>%
  tdmore(prop=0.23) %>%
  estimate(regimen=regimen)

fit <- model %>%
  tdmore(prop=0.23) %>%
  estimate(observed, regimen)

confint(fit)
profile <- profile(fit, maxpts=20)
ggplot(profile, aes(x=ETA1, y=ETA2, z=logLik)) + geom_contour()

newdata = data.frame(TIME=seq(0, 12, length.out=50), CONC=NA)
ggplot(fit %>% predict(newdata), aes(x=TIME, y=CONC)) +
  geom_line(aes(color="Fit")) +
  geom_point(aes(color="Observed"), data=observed) +
  geom_line(aes(color="Population"), data=pred %>% predict(newdata))

ggplot()

stop("STOP")



trueParam <- c(ETA1=0.4, ETA2=-0.3)
true <- myModel %>% predict(times=seq(0, 48), estimate=trueParam)
trueRE <- myModel %>% predict.re.ci(times=seq(0, 48), estimate=trueParam)
observed <- myModel %>% predict.re(estimate=trueParam, times=seq(12, 24)) %>% select(time, CONC)

myEstimate <- myModel %>% estimate(observed)

llSurface <- myModel %>% llSurface(observed=observed)
ggplot(llSurface, aes(x=ETA1, y=ETA2)) + geom_contour(aes(z=V1), bins=400) +
  geom_point(data=trueParam %>% as.list %>% as.data.frame, aes(col="True")) +
  geom_point(x=0, y=0) +
  geom_point(data=myEstimate$res %>% as.list %>% as.data.frame, aes(col="Estimate"))

ipred <- predict.EstimationResult(myEstimate, times=seq(0, 48))
ipredMC <- predictMC.EstimationResult(myEstimate, times=seq(0, 48))

q95 <- function(x) quantile(x, 0.95)
q05 <- function(x) quantile(x, 0.05)
ggplot(true, aes(x=time, y=CONC)) + geom_line(aes(col="True")) +
  geom_line(data=ipred, aes(col="Pred")) +
  geom_point(data=observed) +
  geom_ribbon(data=trueRE, aes(ymin=CONC.CIlower, ymax=CONC.CIupper, fill="True"), alpha=0.5) +
  stat_summary(data=ipredMC, geom="ribbon", fun.ymax=q95, fun.ymin=q05, aes(fill="Pred"), alpha=0.5)
