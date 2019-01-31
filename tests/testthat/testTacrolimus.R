library(tdmore)
library(nlmixr)
library(ggplot2)
library(testthat)
library(dplyr)

#______________________________________________________________________________________________
#----                           TACROLIMUS NLMIXR WITH IOV                                 ----
#______________________________________________________________________________________________

regimen <- data.frame(
 TIME=seq(0, 6)*24,
 AMT=5,
 OCC=c(1,1,2,2,3,3,4)
)

covariates=c(WT=70, HT=1.8, FEMALE=0, CYP3A5=0, PredDose=50, FirstDay=0, HCT=0.45)

plot(tacrolimus_storset, regimen, covariates=covariates, newdata=data.frame(TIME=seq(0, 7*24, by=1), Cwb=NA))

data <- predict(tacrolimus_storset, newdata=seq(0, 7*24, by=1), regimen=regimen, covariates=covariates, parameters = c(EKa_IOV=rnorm(1), EKa_IOV=rnorm(1), EKa_IOV=rnorm(1), EKa_IOV=rnorm(1)))
ggplot(data, aes(x = TIME, y=Cwb)) + geom_line()

observed <- data.frame(TIME=3*24, Cwb=0.007)
ipred <- estimate(tacrolimus_storset, regimen=regimen, covariates=covariates, observed=observed)

## Prediction at day 7 will be quite close to population, but prediction at day 4 should be more accurate
plot(ipred, newdata=data.frame(TIME=seq(0, 7*24, by=1), Cwb=NA))

#______________________________________________________________________________________________
#----                         TACROLIMUS NLMIXR WITHOUT IOV                                ----
#______________________________________________________________________________________________

tacrolimus_nlmixr_no_iov <- nlmixrUI(function(){
  ini({
    TVKa <- 1.01
    #TVTLag <- 0.41  #to be included in model definition itself...

    TVCLp = 811
    TVCLp_CYP3A5 = 1.30
    TVV1p = 6290
    TVQp=1200
    TVV2p=32100

    TVF_PredEMax = 0.67 #relative
    TVF_Pred50 = 35 #mg
    TVF_FirstDay = 2.68 #relative
    TVF_CYP3A5  = 0.82

    # var <- list(CL=0.14842, V1=0.2558818, Q=0.3342555)
    #40%, 54%, 63% resp; 0.43 corrCL/V1, 0.62 corrCL/Q
    ECL + EV1 + EQ ~ c(0.14842,
                       0.08379814, 0.2558818,
                       0.1380948,  0,         0.3342555)
    # Corresponding variance calculated using CV = sqrt( exp(omega^2) -1 )
    # omega^2 = log( CV^2 + 1 )
    EF ~ 0.281337 #57%

    EPS_PROP = 0.149 #standard deviation
  })
  model({
    BMI = WT / HT^2 ;
    FFM = 9.27E3 * WT / (6.68E3 + 216*BMI);
    if(FEMALE) FFM=9.27E3 * WT / (8.78E3 + 244*BMI);
    Ka <- TVKa

    CLp <- TVCLp * (FFM/60)^0.75 * TVCLp_CYP3A5^CYP3A5 * exp(ECL)
    V1p <- TVV1p * exp(EV1)
    V2p <- TVV2p
    Qp <- TVQp * exp(EQ)

    K12 <- Qp/V1p
    K21 <- Qp/V2p
    Ke <- CLp/V1p

    #Covariates: PredDose, FirstDay, CYP3A5
    TVF = 1 * (1 - TVF_PredEMax*PredDose / (PredDose + TVF_Pred50)) * TVF_FirstDay^FirstDay * TVF_CYP3A5^CYP3A5
    F = TVF * exp(EF)

    d/dt(abs) = -Ka*abs
    d/dt(center) = F*Ka*abs - Ke*center - K12*center + K21*periph
    d/dt(periph) = K12*center - K21*periph

    Bmax = 418 #ug/L erythrocytes
    Kd=3.8 #ug/L plasma
    Cp = center/V1p #free concentration in plasma
    Crbc = Cp * HCT * Bmax / (Kd + Cp) #concentration in red blood cells
    Cwb = Cp + Crbc #concentration in whole blood
    Cwb ~ prop(EPS_PROP)
  })
}) %>% tdmore()

regimen <- data.frame(
  TIME=seq(0, 6)*24,
  AMT=5
)

covariates=c(WT=70, HT=1.8, FEMALE=0, CYP3A5=0, PredDose=50, FirstDay=0, HCT=0.45)
plot(tacrolimus_nlmixr_no_iov, regimen, covariates=covariates, newdata=data.frame(TIME=seq(0, 7*24, by=1), Cwb=NA))

data <- predict(tacrolimus_nlmixr_no_iov, newdata=seq(0, 7*24, by=1), regimen=regimen, covariates=covariates)
ggplot(data, aes(x = TIME, y=Cwb)) + geom_line()

observed <- data.frame(TIME=3*24, Cwb=0.007)
ipred <- estimate(tacrolimus_nlmixr_no_iov, regimen=regimen, covariates=covariates, observed=observed)

## Prediction at day 7 will be quite close to population, but prediction at day 4 should be more accurate
plot(ipred, newdata=data.frame(TIME=seq(0, 7*24, by=1), Cwb=NA))

#______________________________________________________________________________________________
#----                           TACROLIMUS ALGEBRAIC WITH IOV                              ----
#______________________________________________________________________________________________

tacrolymusAlgebraicFunction <- function(t, TIME, AMT, ECL, EV1, EQ, EF, EF_IOV, EKa_IOV, WT, HT, FEMALE, CYP3A5, PredDose, FirstDay, HCT) {
  TVKa <- 1.01

  TVCLp = 811
  TVCLp_CYP3A5 = 1.30
  TVV1p = 6290
  TVQp=1200
  TVV2p=32100

  TVF_PredEMax = 0.67 #relative
  TVF_Pred50 = 35 #mg
  TVF_FirstDay = 2.68 #relative
  TVF_CYP3A5  = 0.82

  Bmax = 418 #ug/L erythrocytes
  Kd=3.8 #ug/L plasma

  BMI = WT / HT^2 ;
  FFM = 9.27E3 * WT / (6.68E3 + 216*BMI);
  if(FEMALE) FFM=9.27E3 * WT / (8.78E3 + 244*BMI);
  Ka <- TVKa * exp(EKa_IOV)

  CLp <- TVCLp * (FFM/60)^0.75 * TVCLp_CYP3A5^CYP3A5 * exp(ECL)
  V1p <- TVV1p * exp(EV1)
  V2p <- TVV2p
  Qp <- TVQp * exp(EQ)

  #Covariates: PredDose, FirstDay, CYP3A5
  TVF = 1 * (1 - TVF_PredEMax*PredDose / (PredDose + TVF_Pred50)) * TVF_FirstDay^FirstDay * TVF_CYP3A5^CYP3A5
  F = TVF * exp(EF + EF_IOV)

  K12 <- Qp/V1p
  K21 <- Qp/V2p
  Ke <- CLp/V1p

  Cp <- pk2cptoral_(t, TIME, AMT=F*AMT, V=V1p, KA=Ka, K=Ke, K12=K12, K21=K21)
  Crbc = Cp * HCT * Bmax / (Kd + Cp) #concentration in red blood cells
  return(Cp + Crbc) #concentration in whole blood
}

tacroAlgebraic <- algebraic(tacrolymusAlgebraicFunction)
omega <- matrix(data = c(0.14842,0.08379814,0.1380948,0,0,0,0.08379814,0.2558818,0,0,0,0,0.1380948,0,0.3342555,0,0,0,0,0,0,0.281337,0,0,0,0,0,0,0.05154826,0,0,0,0,0,0,0.891998), nrow = 6, ncol = 6)
colnames(omega) <- c("ECL", "EV1", "EQ", "EF", "EF_IOV", "EKa_IOV")
rownames(omega) <- c("ECL", "EV1", "EQ", "EF", "EF_IOV", "EKa_IOV")

tacroTdm <- tdmore(tacroAlgebraic, omega = omega, res_var = list(errorModel("CONC", prop=0.149)), iov=c("EKa_IOV", "EF_IOV"))

regimen <- data.frame(
  TIME=seq(0, 6)*24,
  AMT=5,
  OCC=seq(1, 7)
)

covariates=c(WT=70, HT=1.8, FEMALE=0, CYP3A5=0, PredDose=50, FirstDay=0, HCT=0.45)

plot(tacroTdm, regimen, covariates=covariates, newdata=data.frame(TIME=seq(0, 7*24, by=1), CONC=NA))

# data <- predict(tacroTdm, newdata=seq(0, 7*24, by=1), regimen=regimen, covariates=covariates)
# ggplot(data, aes(x = TIME, y=CONC)) + geom_line()
#
# data <- predict(tacroTdm, newdata=seq(0, 7*24, by=1), regimen=regimen, covariates=covariates, parameters = c(EKa_IOV=rnorm(1), EKa_IOV=rnorm(1), EKa_IOV=rnorm(1), EKa_IOV=rnorm(1)))
# ggplot(data, aes(x = TIME, y=CONC)) + geom_line()

observed <- data.frame(TIME=c(1,2,3)*24, CONC=c(0.006, 0.0065, 0.007))
ipred <- estimate(tacroTdm, regimen=regimen, covariates=covariates, observed=observed)


DOSE <- findDose(ipred, regimen=regimen, doseRows=c(4), target=data.frame(TIME=4*24, CONC=13.5E-3))
#plot(DOSE)
plot(ipred, DOSE$regimen, newdata=data.frame(TIME=seq(0, 7*24, by=1), CONC=NA)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=12E-3, ymax=15E-3, fill="lightgreen", alpha=0.3)

options(scipen=8)
vcov(ipred)



flatten(data.frame(TIME=10, AMT=10, ADDL=8, II=24))
