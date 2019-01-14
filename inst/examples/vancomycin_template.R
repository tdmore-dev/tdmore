library(tdmore)
library(tidyverse)

## Vancomycin model: 2cpt model with IIV on CL and V2
## Source: Sanchez, J. L., et al. "Population pharmacokinetics of vancomycin in adult and geriatric patients: comparison of eleven approaches." International journal of clinical pharmacology and therapeutics 48.8 (2010): 525-533.
model <- algebraic(function(t, TIME, AMT, ECL, EV2, eGFR, AGE, BW) {
  THETA1=0.157
  THETA5=0.563
  THETA4=0.111
  THETA2=0.283
  THETA3=32.2

  eGFR = eGFR / 1000 * 60 #convert eGFR to L/h
  TVCL = THETA1 + THETA5 * eGFR
  TVV1 = THETA2 * BW
  TVV2 = THETA3 * AGE/53.5
  TVQ = THETA4

  CL = TVCL * exp(ECL) #L/h
  V1 = TVV1 #L
  V2 = TVV2 * exp(EV2) #L
  Q = TVQ #L/h

  # The infusion duration is based on the amount
  DUR <- 1 #1 hour by default
  if(AMT >= 1000) DUR <- 2 #high doses are infused longer

  pk2cptivinfusion()
}) %>% tdmore(omega=c(ECL=0.2449^2, EV2=0.068^2),
              res_var=errorModel(prop=0.249))

## Our targeted PK metric
TARGET <- 15 #in this case, we target Ctrough [mg/L]
SAFETY <- 50 #limit for peak post-infusion  [mg/L]

## We can now repeat the exercise for a different patient
covariates <- c(eGFR=125, AGE=45, BW=70)   ### YOU NEED TO ADAPT THIS
LOADING_DOSE <- data.frame(TIME=0, II=0, ADDL=0, AMT=25*covariates['BW'])
regimen <- bind_rows(
  LOADING_DOSE,
  data.frame(TIME=12, II=12, ADDL=0, AMT=1000) ### YOU NEED TO ADAPT THIS
)
observed <- bind_rows(
  data.frame(TIME=24, CONC=7.6) ### YOU NEED TO ADAPT THIS
)

## A priori plot
plot(model, regimen=regimen, covariates=covariates) +
  geom_point(aes(y=CONC), data=observed) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  scale_x_continuous(breaks=seq(0, 10*24, by=12))

## Estimate
ipred <- estimate(model, regimen=regimen, covariates=covariates, observed=observed)
plot(ipred) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  scale_x_continuous(breaks=seq(0, 10*24, by=12))

## Predict the dose
regimen <- bind_rows(
  LOADING_DOSE,
  data.frame(TIME=12, II=12, ADDL=0, AMT=1000), #the past   ### YOU NEED TO ADAPT THIS
  data.frame(TIME=24, II=12, ADDL=11, AMT=NA) # the future
)
recommendation <- findDose(ipred, regimen, target=data.frame(TIME=4*24, CONC=TARGET))
print(recommendation)
plot(ipred, regimen=recommendation$regimen, population=F) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  scale_x_continuous(breaks=seq(0, 10*24, by=12))

