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
if(exists("myRxODEModel")) myRxODEModel$dynUnload()
if(exists("myRxODEModel") && myRxODEModel$model == modelCode) { #skip
} else {
  myRxODEModel <- RxODE(modelCode)
}

myModel <- Model(
  predict=function(times, estimates=c(ETA1=0, ETA2=0)) {
    observed <- data.frame(time=times)
    ev <- eventTable()
    ev$add.dosing(dose=5, nbr.doses=10, dosing.interval = 24)
    ev$add.sampling(time=observed$time)
    res <- myRxODEModel$solve(params=c(observed[1,], estimates), events=ev, inits=c(0, 0, 0))
    if(!is.matrix(res)) {
      res <- matrix(res, nrow=1, dimnames = list(NULL, names(res)))
    }
    res %>% as.data.frame()
  },
  parameters=myRxODEModel$get.modelVars()$params,
  propSigma = 0.23 #23%
)

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
