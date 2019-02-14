##
## Purpose of script:
## Demonstrate model-based closed loop control
##
## Author: Ruben Faelens
##
## Date Created: Wed Feb 06 15:51:42 2019
##
## Copyright (c) Ruben Faelens, 2019
## Email: ruben.faelens@gmail.com
##
## ---------------------------
##
## Notes:
## Closed-loop control is a technique where the fit is adapted from day-to-day, each time
## a new concentration comes in.
##
## It does this by assuming the past is fixed.
## Furthermore, the current EBE priors are based on the previous EBE estimate, and its uncertainty.
##
## This is implemented in TDMore by doing the following:
## - Define a model with only IIV. Mark the IIV as IOV.
## - The model should take the typical values from the covariates.
## 1) Estimate EBE for the first occasion. Use the population typical values and population OMEGA.
## regimen <- regimen %>% filter(OCC==1)
## theta <- c(CL=18, V1=28, V2=309, Q=20) #true THETA
## estimate(m1, observed=observed, regimen=regimen, covariates=theta)
## 2) For the next occasion, estimate as follows:
## m1$omega <- vcov(previousIpred)
## theta <- predict(previousIpred, newdata=0)[, c("CL", "V1", "V2", "Q")]
## names(theta) <- c("TVCL", "TVV1", "TVV2", "TVQ")
## fixedParameters <- coef(previousIpred)
## estimate(m1, observed=observed, regimen=regimen, covariates=theta, par=fixedParameters)
## 3) For the next occasion, estimate as follows:
## m1$omega <- vcov(previousIpred)
## theta <- predict(previousIpred, newdata=24)[, c("CL", "V1", "V2", "Q")]
## names(theta) <- c("TVCL", "TVV1", "TVV2", "TVQ")
## fixedParameters <- c(fixedParameters, coef(previousIpred))
## estimate(m1, observed=observed, regimen=regimen, covariates=theta, par=fixedParameters)
##
## ---------------------------
library(nlmixr)
library(tdmore)
library(tidyverse)
## Special adaptations in the model
## 1) Pass THETA values as covariates
## 2) Mark all IIV as IOV
## 3) Return parameter values in output
m1 <- algebraic(function(t, TIME, AMT, EKa, ECL, EV1, EQ, EV2, TVKa, TVCL, TVV1, TVQ, TVV2){
    Ka <- TVKa*exp(EKa)
    CL <- TVCL * exp(ECL)
    V1 <- TVV1 * exp(EV1)
    V2 <- TVV2 * exp(EV2)
    Q <- TVQ * exp(EQ)

    CONC = pk2cptoral_(t, TIME, AMT, V=V1, K=CL/V1, KA=Ka, K12=Q/V1, K21=Q/V2)
    list(
      Cwb = CONC * 1000, #concentration in whole blood, in ng/mL or ug/L
      Ka=Ka, CL=CL, V1=V1, V2=V2, Q=Q
      )
}, output="Cwb") %>% tdmore(
  omega=c(
    EKa=0.22,
    ECL=0.44,
    EV1=0.62,
    EQ=0.36,
    EV2=0.84
  )**2,
  iov=c("EKa", "ECL", "EV1", "EQ", "EV2"),
  res_var=errorModel(var="Cwb", add=0.1, prop=0.15)
)
theta <- data.frame(TIME=0, TVKa=0.5, TVCL=16.1, TVV1=133, TVQ=68.3, TVV2=348)

# db <- readRDS("D:/OtterDrive/Tacrolimus/3.Analysis Data/PostTransplant_Vanhove2017BJCP_KWS.Rds")
# subject <- db %>% ungroup() %>% filter(ID==1) %>% filter(
#   (EVID == 0 & MDV.all == 0) | EVID==1
# ) %>% select(-(Age:TxDate), -FORM, -MDV.censor, -MDV.all, -MDV, -TxNR, -ID, -DAY, -TATx, -TSLD)
# subject %>% mutate_all(as.numeric) %>% as.data.frame %>% write.csv(file="", row.names=F)

subject <- '
"OCC","TSFD","CONC","EVID","AMT"
1,0,NA,1,10
1,22.4,9,0,NA
2,22.45,NA,1,5
2,34.45,NA,1,6
2,46.4,17,0,NA
3,46.45,NA,1,6
3,58.45,NA,1,5.5
3,70.4,21,0,NA
4,70.45,NA,1,5.5
4,82.45,NA,1,0
4,94.4,13,0,NA
5,94.45,NA,1,3
5,106.45,NA,1,2.5
5,118.4,10,0,NA
6,118.45,NA,1,3
6,130.45,NA,1,3.5
6,142.4,9,0,NA
7,142.45,NA,1,3.5
7,154.45,NA,1,4.5
7,166.4,7,0,NA
8,166.45,NA,1,4
8,178.45,NA,1,6
8,190.4,9,0,NA
9,190.45,NA,1,5
9,202.45,NA,1,6
9,214.4,13,0,NA
10,214.45,NA,1,5.5
10,226.45,NA,1,5.5
10,238.4,14,0,NA
11,238.45,NA,1,5.5
11,250.45,NA,1,5.5
11,262.4,15,0,NA
12,262.45,NA,1,11
12,286.4,10,0,NA
13,286.45,NA,1,11
13,310.4,10,0,NA
14,310.45,NA,1,11
14,334.4,10,0,NA
15,334.45,NA,1,12
' %>% read_csv

regimen <- subject %>% filter(EVID==1) %>%
  transmute(TIME=TSFD, AMT, OCC) %>% as.data.frame
observed <- subject %>% filter(EVID==0) %>% transmute(
  TIME=TSFD-0.05, #make sure it is before the administration...
  Cwb=CONC
) %>% as.data.frame
covariates <- theta ##TODO: time-varying covariates

## The truth is probably not right; you allow the IIV to change every day!
#truth <- estimate(m1, regimen=regimen, observed=observed, covariates=covariates, se.fit=FALSE,
#                control=list(trace=1)) # estimate full profile

N <- length(m1$iov) #how many IOV parameters there are

## First iteration
i <- 1
fixedParameters <- c()
covariates <- theta
input <- observed[i,]
#debugonce(tdmore:::plot.tdmorefit)
ipred <- estimate(m1, regimen=regimen %>% filter(OCC <= 1), observed=input, covariates=covariates, # Add filter on OCC
                  control=list(trace=1, REPORT=1))
plot(ipred)

## Second iteration
for(i in 2:5) {

  previousEta <- coef(ipred)[  seq(1, N) + (i-2)*N ] #fix the previous Eta
  fixedParameters <- c(fixedParameters, previousEta)

  previousEbe <- predict(ipred, newdata=observed$TIME[i-1])[, c("Ka", "CL", "V1", "V2", "Q")] #set up a new THETA
  names(previousEbe) <- paste0("TV", names(previousEbe)) #use these as the typical values for the new estimation
  previousEbe$TIME <- observed$TIME[i-1]

  ## TODO: take the right rows from vcov (the ones that refer to the uncertainty on the current occasion)
  m1$omega <- vcov(ipred)[seq(1,N)+(i-2)*N, seq(1,N)+(i-2)*N] #adapt OMEGA

  input <- observed[i,]
  covariates <- bind_rows(covariates, previousEbe)
  # So in this call, we use the estimated ETA and real THETA for occasion 1
  # And in occasion 2, we estimate the new ETA using an OMEGA of vcov, and the THETA that is the previous EBE

  #debugonce(tdmore:::estimate)
  ipred <- estimate(m1, regimen=regimen %>% filter(OCC <= i), observed=input, covariates=covariates,
                    fix=fixedParameters,  ## TODO: add a mechanism to FIX parameters
                    control=list(trace=1, REPORT=1))
  #debugonce(tdmore:::plot.tdmorefit)
  print(plot(ipred, se.fit=NA))
}
