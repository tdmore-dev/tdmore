library(tdmore)
library(testthat)
library(RxODE)
library(magrittr)

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

regimen <- data.frame(
  TIME=seq(0, 3)*24,
  AMT=5 #5mg
)

expect_snapshot_output(
  print(tdmore)
)

expect_snapshot_output(
  print(summary(tdmore))
)

tdmore:::model_predict.RxODE(m1, times=seq(0,16), regimen=regimen, parameters=c(ECL=0, EVc=0))

# Default tdmore plot
tdmore:::predict.tdmore(tdmore, newdata=seq(0,16), regimen=regimen)
plot(tdmore, regimen=regimen, se.fit=F)

# Compute PRED
pred <- tdmore %>% tdmore:::estimate(regimen = regimen)
expect_equal(pred$res, c(EVc=0.0, ECL=0.0))

# Compute IPRED
observed <- data.frame(TIME=c(2), CONC=c(0.040))
ipred <- tdmore %>% tdmore:::estimate(observed = observed, regimen = regimen)
expect_equal(round(ipred$res, digits=4), c(EVc=0.1175, ECL=0.0336))

# Default IPRED plot
plot(ipred, se.fit=F) + ggplot2::coord_cartesian(xlim=c(0, 100))

# Test combinations of ADDL and SS and II
expect_snapshot_value(
  predict(as.population(tdmore), newdata=seq(0, 14)),
  "serialize"
)
expect_snapshot_value(
  predict(as.population(tdmore), newdata=seq(0, 14), regimen=data.frame(TIME=0, AMT=100)),
  "serialize"
)

## RxODE throws an error 'ii requires non zero additional doses'
## but tdmore does support this treatment regimen,
## as it is useful to plan e.g. additional doses
predict(as.population(tdmore), newdata=seq(0, 14), regimen=data.frame(TIME=0, AMT=100, II=12))

## eventTable$add.regimen would give an error
## "ii requires non zero additional doses"
## But a manual event table works perfectly!
pop <- as.population(tdmore)
predict <- function(..., digits=1e-4) {
  result <- stats::predict(...)
  round(result, digits=digits)
}
expect_snapshot_value(
  predict(pop, newdata=seq(0, 14), regimen=data.frame(TIME=0, AMT=100, II=12, ADDL=1)),
  "serialize"
)

tdmoreMonolixSS <- m1 %>%
  tdmore(omega=omegas,
         res_var=list(errorModel("CONC", prop=0.23)), #Model has 23% proportional error
         nbSSDoses = 5
  )
expect_snapshot_value(
  predict(as.population(tdmoreMonolixSS), newdata=seq(0, 14)+12*10, regimen=data.frame(TIME=0+12*10, AMT=100, II=12, SS=1)),
  "serialize"
)

expect_snapshot_value(
  predict(pop, newdata=seq(0, 14), regimen=data.frame(TIME=0, AMT=100, II=12, SS=1)),
  "serialize",

)

expect_snapshot_value(
  predict(pop, newdata=seq(0, 14), regimen=data.frame(TIME=0, AMT=100, II=12, SS=1, ADDL=0)),
  "serialize"
)

expect_error(
 predict(pop, newdata=seq(0, 14), regimen=data.frame(TIME=0, AMT=100, II=12, SS=1, ADDL=1)),
 regexp = "'ss' with 'addl' not supported"
)


# Test IOV with two treatments at same time -------------------------------
time=15
regimen=data.frame(TIME=c(0,0, 24, 24, 48), AMT=5, OCC=seq(1, 5))
covariates=NULL

m1ModelIOV <- m1 %>%
  tdmore(omega=omegas,
         res_var=list(errorModel("CONC", prop=0.23)),
         iov=names(omegas)) #Model has 23% proportional error

## This should work fine!
predict(as.population(m1ModelIOV, regimen=regimen), newdata=time)

