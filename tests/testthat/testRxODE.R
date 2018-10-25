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

tdmore <- RxODE::RxODE(modelCode) %>%
  tdmore(omega=vectorToDiagonalMatrix(omegas), prop=0.23) #Model has 23% proportional error

regimen <- data.frame(
  TIME=seq(0, 3)*24,
  AMT=5 #5mg
)

# Default tdmore plot
plot(tdmore, regimen, vars=c("perip", "centr"))

# Compute PRED
pred <- tdmore %>% estimate(regimen = regimen)
expect_equal(pred$res, c(ECL=0.0, EVc=0.0))

# Compute IPRED
observed <- data.frame(TIME=c(2), CONC=c(0.040))
ipred <- tdmore %>% estimate(observed = observed, regimen = regimen)
expect_equal(round(ipred$res, digits=4), c(ECL=0.0336, EVc=0.1175))

# Default IPRED plot
plot(ipred)

