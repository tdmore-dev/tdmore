library(tdmore)
library(testthat)
library(RxODE)
library(magrittr)

context("Test that the RxODE model class works as intended")

modelCode <- "
CL = 3.7 * exp(ECL);
Vc = 61 * exp(EVc);
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
omegas=c(EVc=0.19^2, ECL=0.28^2)

model <- RxODE::RxODE(modelCode) %>%
  tdmore(omega=convertVectorToDiag(omegas), prop=0.23) #Model has 23% proportional error

regimen <- data.frame(
  TIME=seq(0, 7)*24,
  AMT=5 #5mg
)

observed <- data.frame(TIME=c(2), CONC=c(0.040))

# Compute PRED
pred <- model %>% estimate(regimen = regimen)
expect_equal(pred$res, c(ECL=0.0, EVc=0.0))

# Compute IPRED
ipred <- model %>% estimate(observed = observed, regimen = regimen)
expect_equal(round(ipred$res, digits=4), c(ECL=0.0336, EVc=0.1175))

# Custom plot, compare PRED & IPRED
library(ggplot2)
ggplot(observed, aes(x=TIME, y=CONC)) + geom_point() +
  geom_line(aes(color="Population"), data=predict(pred, newdata=seq(0, 48, by=0.1))) +
  geom_line(aes(color="Individual"), data=predict(ipred, newdata=seq(0, 48, by=0.1)))

# Default TDMore plot for IPRED
plot(ipred, newdata=data.frame(TIME=seq(0, 48, by=0.1), CONC=NA))
