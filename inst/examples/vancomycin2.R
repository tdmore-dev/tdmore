library(tdmore)
library(tidyverse)

#______________________________________________________________________________________________
#----                                   VANCOMYCIN MODEL                                   ----
#______________________________________________________________________________________________
# Vancomycin model: 2cpt model with IIV on CL and V2
# Source: Sanchez, J. L., et al. "Population pharmacokinetics of vancomycin in adult and geriatric patients: comparison of eleven approaches."
# International journal of clinical pharmacology and therapeutics 48.8 (2010): 525-533.

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

#______________________________________________________________________________________________
#----               EXERCICE 1a: Plotting regimen               ----
#______________________________________________________________________________________________

# Our targeted PK metric
TARGET <- 15 # in this case, we target Ctrough [mg/L]
SAFETY <- 50 # limit for peak post-infusion  [mg/L]
STEADY_STATE <- 72 # steady state is reached after 3 days

# Patient covariates
covariates <- c(eGFR=60, AGE=73, BW=55)
LOADING_DOSE <- data.frame(TIME=0, II=0, ADDL=0, AMT=25*covariates['BW'] )

regimen <- bind_rows(
  LOADING_DOSE,                                # Past
  data.frame(TIME=12, II=0, ADDL=0, AMT=750),  # Past
  data.frame(TIME=24, II=12, ADDL=4, AMT=1000) # Past (1st recommendation)
)

observed <- data.frame(TIME=c(24, 72), CONC=c(12.5, 17.5))

plot(model, regimen=regimen, covariates=covariates) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  geom_vline(xintercept=72, linetype=2) +
  geom_point(aes(x=TIME, y=CONC), data=observed) +
  scale_x_continuous(breaks=seq(0, 84, by=12))

#______________________________________________________________________________________________
#---- EXERCICE 1b: Optimisation based on first two observations and forward predictions    ----
#______________________________________________________________________________________________

ipred <- estimate(model, observed = observed, covariates = covariates, regimen = regimen)

recommendation <- findDose(ipred, regimen=regimen, target=data.frame(TIME=STEADY_STATE, CONC=TARGET))
print(recommendation)

observed <- data.frame(TIME=c(24, 72, 6*24), CONC=c(12.5, 17.5, 20.2))

regimen <- bind_rows(
  LOADING_DOSE,                                # Past
  data.frame(TIME=12, II=0, ADDL=0, AMT=750),  # Past
  data.frame(TIME=24, II=12, ADDL=10, AMT=1000) # Past (1st recommendation)
)

plot(model, regimen=regimen, covariates=covariates) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  geom_vline(xintercept=144, linetype=2) +
  geom_point(aes(x=TIME, y=CONC), data=observed) +
  annotate("text", x=144, y=TARGET, label=TARGET, vjust=-0.5, hjust=1.5) +
  annotate("text", x=144, y=observed$CONC[3], label=observed$CONC[3], vjust=-0.5, hjust=1.5)
  #scale_x_continuous(breaks=seq(0, 84, by=12))

#########

observed <- data.frame(TIME=c(6*24), CONC=c(20.2))

regimen <- bind_rows(
  LOADING_DOSE,                                   # Past
  data.frame(TIME=12, II=0, ADDL=0, AMT=750),     # Past
  data.frame(TIME=24, II=12, ADDL=9, AMT=1000)   # Past (1st recommendation)
)

ipred <- estimate(model, observed = observed, covariates = covariates, regimen = regimen)

regimen <- bind_rows(
  LOADING_DOSE,                                   # Past
  data.frame(TIME=12, II=0, ADDL=0, AMT=750),     # Past
  data.frame(TIME=24, II=12, ADDL=9, AMT=1000),   # Past (1st recommendation)
  data.frame(TIME=6*24, II=12, ADDL=14, AMT=NA)  # Past (2st recommendation)
)

recommendation <- findDose(ipred, regimen=regimen, target=data.frame(TIME=6*24+STEADY_STATE, CONC=TARGET))
print(recommendation)

OBSERVED <- data.frame(TIME=c(6*24, 7*24), CONC=c(20.2, 17.1))
plot(ipred, regimen=recommendation$regimen, covariates=covariates) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  geom_vline(xintercept=168, linetype=2) +
  scale_x_continuous(breaks=seq(0, 240, by=12), limits = c(0, 240)) +
  scale_y_continuous(limits = c(5, 80)) +
  geom_point(aes(x=TIME, y=CONC), data = OBSERVED, size=I(3))

#########################

regimen <- bind_rows(
  LOADING_DOSE,                                   # Past
  data.frame(TIME=12, II=0, ADDL=0, AMT=750),     # Past
  data.frame(TIME=24, II=12, ADDL=9, AMT=1000),   # Past (1st recommendation)
  data.frame(TIME=6*24, II=8, ADDL=14, AMT=NA)  # Past (2st recommendation)
)

recommendation <- findDose(ipred, regimen=regimen, target=data.frame(TIME=6*24+STEADY_STATE, CONC=TARGET))
print(recommendation)

# plot(ipred, regimen=recommendation$regimen, covariates=covariates) +
#   geom_hline(yintercept=c(TARGET, SAFETY)) +
#   geom_vline(xintercept=0, linetype=2) +
#   scale_x_continuous(breaks=seq(0, 240, by=12), limits = c(0, 240)) +
#   scale_y_continuous(limits = c(5, 80))

#########################


OBSERVED <- data.frame(TIME=c(6*24, 7*24), CONC=c(20.2, 17.1))

plot(ipred, regimen=recommendation$regimen, covariates=covariates) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  geom_vline(xintercept=168, linetype=2) +
  scale_x_continuous(breaks=seq(0, 240, by=12), limits = c(0, 240)) +
  scale_y_continuous(limits = c(5, 80)) +
  geom_point(aes(x=TIME, y=CONC), data = OBSERVED, size=I(3))
