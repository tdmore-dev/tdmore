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
m1 <- RxODE::RxODE(modelCode)
tdmore <- m1 %>%
  tdmore(omega=omegas,
         res_var=list(errorModel("CONC", prop=0.23))) #Model has 23% proportional error

tdmoreMonolixSS <- m1 %>%
  tdmore(omega=omegas,
         res_var=list(errorModel("CONC", prop=0.23)), #Model has 23% proportional error
          nbSSDoses = 5
  )

regimen <- data.frame(
  TIME=seq(0, 3)*24,
  AMT=5 #5mg
)

expect_known_output(
  print(tdmore),
  "tdmoreRxode.txt"
)

expect_known_output(
  print(summary(tdmore)),
  "tdmoreRxodeSummary.txt"
)

model_predict(m1, times=seq(0:16), regimen=regimen, parameters=c(ECL=0, EVc=0))

# Default tdmore plot
plot(tdmore, regimen)

# Compute PRED
pred <- tdmore %>% estimate(regimen = regimen)
expect_equal(pred$res, c(ECL=0.0, EVc=0.0))

# Compute IPRED
observed <- data.frame(TIME=c(2), CONC=c(0.040))
ipred <- tdmore %>% estimate(observed = observed, regimen = regimen)
expect_equal(round(ipred$res, digits=4), c(ECL=0.0336, EVc=0.1175))

# Default IPRED plot
plot(ipred)

# Test combinations of ADDL and SS and II
expect_known_value(
  predict(tdmore, newdata=seq(0, 14)),
  "Rxode.empty.regimen.prediction"
)
expect_known_value(
  predict(tdmore, newdata=seq(0, 14), regimen=data.frame(TIME=0, AMT=100)),
  "Rxode.regimen.1dose"
)

## RxODE throws an error 'ii requires non zero additional doses'
## but tdmore does support this treatment regimen,
## as it is useful to plan e.g. additional doses
predict(tdmore, newdata=seq(0, 14), regimen=data.frame(TIME=0, AMT=100, II=12))

expect_error(
  predict(tdmore, newdata=seq(0, 14), regimen=data.frame(TIME=0, AMT=100, II=12, ADDL=0)),
  regexp="ii requires non zero additional doses"
)
expect_known_value(
  predict(tdmore, newdata=seq(0, 14), regimen=data.frame(TIME=0, AMT=100, II=12, ADDL=1)),
  "Rxode.regimen.2dose"
)

expect_known_value(
  predict(tdmoreMonolixSS, newdata=seq(0, 14)+12*10, regimen=data.frame(TIME=0+12*10, AMT=100, II=12, SS=1)),
  "Rxode.regimen.1doseSSMonolixStyle"
)

expect_known_value(
  predict(tdmore, newdata=seq(0, 14), regimen=data.frame(TIME=0, AMT=100, II=12, SS=1)),
  "Rxode.regimen.1doseSS"
)

expect_known_value(
  predict(tdmore, newdata=seq(0, 14), regimen=data.frame(TIME=0, AMT=100, II=12, SS=1, ADDL=0)),
  "Rxode.regimen.1doseSS"
)
expect_error(
 predict(tdmore, newdata=seq(0, 14), regimen=data.frame(TIME=0, AMT=100, II=12, SS=1, ADDL=1)),
 regexp = "ss with addl not supported yet"
)


# Test IOV with two treatments at same time -------------------------------
time=15
regimen=data.frame(TIME=c(0,0, 24, 24, 48), AMT=5, OCC=seq(1, 5))
covariates=NULL

m1ModelIOV <- m1 %>%
  tdmore(omega=omegas,
         res_var=list(errorModel("CONC", prop=0.23)),
         iov=names(omegas)) #Model has 23% proportional error

predict(m1ModelIOV, regimen=regimen, newdata=time)

