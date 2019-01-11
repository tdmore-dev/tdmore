library(tdmore)
library(tidyverse)

## Vancomycin model: 2cpt model with IIV on CL and V2
## Source: Sanchez, J. L., et al. "Population pharmacokinetics of vancomycin in adult and geriatric patients: comparison of eleven approaches." International journal of clinical pharmacology and therapeutics 48.8 (2010): 525-533.
model <- algebraic(function(t, TIME, AMT, ECL, EV2, CRCL, AGE, BW) {
  THETA1=0.157
  THETA5=0.563
  THETA4=0.111
  THETA2=0.283
  THETA3=32.2

  CRCL = CRCL / 1000 * 60 #convert CRCL to L/h
  TVCL = THETA1 + THETA5 * CRCL
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

## Patient covariates
covariates <- c(CRCL=125, AGE=45, BW=70)
LOADING_DOSE <- data.frame(TIME=0, II=0, ADDL=0, AMT=25*covariates['BW'] )

## What should the regimen be for this patient, based on the dosing table?
regimen_dosingTable <- bind_rows(
  LOADING_DOSE,
  data.frame(TIME=12, II=12, ADDL=4, AMT=1500)
)

# How does this look like?
# Predict the typical individual (open loop dosing based on patient covariates)
plot(model, regimen=regimen_dosingTable, covariates=covariates) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  geom_vline(xintercept=0, linetype=2) +
  scale_x_continuous(breaks=seq(0, 200, by=12))
  # This is a concentration-time plot
  # It shows (in blue) the estimated value and 95% prediction interval
  # for the covariate-corrected population.

## What was actually administered?
# Dosing regimen: 1000 mg q12
regimen <- bind_rows(
  LOADING_DOSE,
  data.frame(TIME=12, II=12, ADDL=4, AMT=1000)
)

# Predict the typical individual (open loop dosing based on patient covariates)
plot(model, regimen=regimen, covariates=covariates) +
  geom_hline(yintercept=c(TARGET, SAFETY))  +
  geom_vline(xintercept=0, linetype=2) +
  scale_x_continuous(breaks=seq(0, 200, by=12)) +
  annotate("text", x=24, y=TARGET, label=TARGET, vjust=-0.5, hjust=1.5) +
  geom_text( aes(y=CONC, label=round(CONC, 2)), data=predict(model, regimen=regimen, covariates=covariates, newdata=24))
  # you are not hitting your target
  #   is this a clinical issue?
  #   would you like to change the dose?

## Fortunately, we measured concentrations
# We now input the observed concentration from Day 1
observed <- bind_rows(
  data.frame(TIME=24, CONC=7.6)
)
plot(model, regimen=regimen, covariates=covariates) +
  geom_point(aes(x=TIME, y=CONC), data=observed) +
  geom_text(aes(x=TIME, y=CONC, label=CONC), data=observed, hjust=1.5) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  geom_vline(xintercept=24, linetype=2) +
  scale_x_continuous(breaks=seq(0, 200, by=12))
    # the point is in the blue area, so this is a subject that is within the 95% PI for the population
    # If it would not be in the blue band:
    #   perhaps you made a mistake in the input data?
    #   or maybe the model is not appropriate?
    #   or maybe that patient is simply very unusual? (OMEGA usually has large RSE)
    #


# Estimate the patient PK parameters, and get a better prediction
ipred <- estimate(model, regimen=regimen, covariates=covariates, observed=observed)
coef(ipred)

TVCL = 0.157 + 0.563 * covariates['CRCL'] / 1000 * 60
message("Typical clearance is ", TVCL)
CL = TVCL * exp(coef(ipred)['ECL'])
message("Patient estimated clearance is ", CL)

# We can plot this estimate of the patient PK
plot(ipred, se.fit=F) +
  geom_hline(yintercept=c(TARGET, SAFETY))  +
  geom_vline(xintercept=24, linetype=2) +
  scale_x_continuous(breaks=seq(0, 200, by=12))
  # This is a time-concentration plot
  # It shows (in red) the estimated value and 95% prediction interval for the given patient
  # This prediction is more accurate now, because we have a concentration point to characterize the patient.
  #
  # The curve does not go exactly through the point.
  # The bayesian approach penalizes moving estimated parameters away from the typical value.
  # It only allows doing so to fit the given observations.
  #
  # It also shows (in blue) the population typical value and 95% population prediction interval
  # This is the unadjusted prediction for the population (without concentration samples)
  #
  #   Does this estimate describe what you observed?
  #


# Now we can adapt the dose, and predict into the future
# We can evaluate a dose of 2000mg q12 manually
regimen <- bind_rows(
  LOADING_DOSE,
  data.frame(TIME=12, II=12, ADDL=0, AMT=1000), #the past
  data.frame(TIME=24, II=12, ADDL=10, AMT=2000) #the future
)
plot(ipred, regimen=regimen, se.fit=T, population=F) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  geom_vline(xintercept=24, linetype=2) +
  scale_x_continuous(breaks=seq(0, 200, by=12))

# or 1000mg q8
regimen <- bind_rows(
  LOADING_DOSE,
  data.frame(TIME=12, II=12, ADDL=0, AMT=1000), #the past
  data.frame(TIME=24, II=8, ADDL=10, AMT=1000)  #the future
)
plot(ipred, regimen=regimen, se.fit=T, population=F) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  geom_vline(xintercept=24, linetype=2) +
  scale_x_continuous(breaks=seq(0, 200, by=12))


## We can also ask TDMore which dose should be given (for a pre-specified dosing interval)
regimen <- bind_rows(
  LOADING_DOSE,
  data.frame(TIME=12, II=12, ADDL=0, AMT=1000), #the past
  data.frame(TIME=24, II=12, ADDL=10, AMT=NA) #TDMore will always try to adapt the last row in the regimen
)
recommendation <- findDose(ipred, regimen=regimen,
                           target=data.frame(TIME=3*24, CONC=TARGET) #what should we aim for?
                           )
recommendation

plot(ipred, regimen=recommendation$regimen, se.fit=T, population=F) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  geom_vline(xintercept=24, linetype=2) +
  scale_x_continuous(breaks=seq(0, 200, by=12))
  ## The dose we find should be evaluated on safety as well
    # The predicted concentration is too high
    # We should probably switch to a q8 regimen

## We can also ask TDMore which dose should be given every 8h
regimen <- bind_rows(
  LOADING_DOSE,
  data.frame(TIME=12, II=12, ADDL=0, AMT=1000), #the past
  data.frame(TIME=24, II=8, ADDL=10, AMT=NA) #TDMore will always try to adapt the last row in the regimen
)
recommendation <- findDose(ipred, regimen=regimen,
                           target=data.frame(TIME=3*24, CONC=TARGET) #what should we aim for?
)
print(recommendation)
plot(ipred, regimen=recommendation$regimen, se.fit=T, population=F) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  geom_vline(xintercept=24, linetype=2) +
  scale_x_continuous(breaks=seq(0, 200, by=12))

# We can ask which loading dose should be given, and which maintenance dose
regimen <- bind_rows(
  LOADING_DOSE,
  data.frame(TIME=12, II=12, ADDL=0, AMT=1000),
  data.frame(TIME=24, II=8, ADDL=0, AMT=NA) #single loading dose
)
loadingRecommendation <- findDose(ipred, regimen=regimen,
                           target=data.frame(TIME=24+8, CONC=TARGET) )
print(loadingRecommendation)
plot(ipred, regimen=loadingRecommendation$regimen, population=F) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  geom_vline(xintercept=24, linetype=2) +
  scale_x_continuous(breaks=seq(0, 200, by=12))
# and now the maintenance dose
regimen <- bind_rows(
  loadingRecommendation$regimen, #the previously found regimen
  data.frame(TIME=24+8, II=8, ADDL=10, AMT=NA) #maintenance dose
)
print(regimen)
recommendation <- findDose(ipred, regimen=regimen,
                           target=data.frame(TIME=3*24, CONC=TARGET) )
print(recommendation)
print(recommendation$regimen)
plot(ipred, regimen=recommendation$regimen, population=F) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  geom_vline(xintercept=24, linetype=2) +
  scale_x_continuous(breaks=seq(0, 200, by=12))
    ## The loading dose and maintenance dose are so similar
    ## that a loading dose is probably not required here

## How does the TDMore dose compare to the dose chosen by the physician in this patient?
print(recommendation)
    # It is exactly the same

## Were we right with our prediction?
regimen <- bind_rows(
  LOADING_DOSE,
  data.frame(TIME=12, II=12, ADDL=0, AMT=1000),
  #This is the regimen that was chosen by the physician (and TDMore)
  data.frame(TIME=24, II=8, ADDL=5, AMT=1000)
)
predictedValue <- predict(ipred, regimen=regimen, newdata=3*24, se.fit=T)
plot(ipred, regimen=regimen, population=F) +
  annotate("point", x=3*24, y=14.2, shape=2) +
  annotate("text", x=3*24, y=14.2, label=14.2, hjust=1.5) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  geom_vline(xintercept=24, linetype=2) +
  scale_x_continuous(breaks=seq(0, 200, by=12)) +
  geom_text(aes(x=TIME, y=CONC, label=round(CONC,1)),
    data=predictedValue, hjust=1.5 )
## And the actual numbers?
print(predictedValue)
# What was actually observed: 14.2 mg/L


## Limitations
# No time-varying covariates in algebraic equations.
# You would need an ODE model for that (but it works in TDMore)
#
# You cannot search for the most appropriate dosing interval yet.
#
# More advanced 'dose searching' methods are planned.
# Eg find the most appropriate dose + dosing interval to
#   1) exceed the target for Ctrough
#   2) while keeping CMax below a target value
#
# No automated detection of "trough", "steady-state trough", "3-hours post-infusion", etc.
#   Our aim is to allow the following to work:
#   findDose(ipred, regimen=regimen, target=data.frame(TIME=next_trough(), CONC=TARGET) )
