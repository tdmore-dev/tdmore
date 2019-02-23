#### TODO
#### Algebraic equations with IOV will never match with the ODE equations.
#### How can we adapt the algebraic equations to be more correct?
#### TODO


library(tdmore)
library(nlmixr)
library(ggplot2)
library(testthat)
library(dplyr)

regimen <- data.frame(
  TIME=seq(0, 6)*24,
  AMT=5,
  OCC=seq(1, 7)
)

covariates=c(WT=70, HT=1.8, FEMALE=0, CYP3A5=0, PredDose=50, FirstDay=0, HCT=0.45)

observed <- data.frame(TIME=c(1,2,3)*24, Cwb=c(6, 6.5, 7)) # Cwb (see tacrolimus model)

#______________________________________________________________________________________________
#----                            TACROLIMUS (RXODE MODEL)                                 ----
#______________________________________________________________________________________________

tacro_ode <- tacrolimus_storset # Loaded from the model library

# Population plot
plot(tacro_ode, regimen, covariates=covariates, newdata=data.frame(TIME=seq(0, 7*24, by=1), Cwb=NA))

# Compute ipred based on observed data
ipred_ode <- estimate(tacro_ode, regimen=regimen, covariates=covariates, observed=observed)

# Find next dose
dose_ode <- findDose(ipred_ode, regimen=regimen, doseRows=c(4), target=data.frame(TIME=4*24, Cwb=13.5))
plot(ipred_ode, dose_ode$regimen, newdata=data.frame(TIME=seq(0, 7*24, by=1), Cwb=NA), se.fit=T) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=12, ymax=15, fill="lightgreen", alpha=0.3)

#______________________________________________________________________________________________
#----                            TACROLIMUS (ALGEBRAIC MODEL)                              ----
#______________________________________________________________________________________________

tacrolimusAlgebraicFunction <- function(t, TIME, AMT, ECL, EV1, EQ, EF, EF_IOV, EKa_IOV, WT, HT, FEMALE, CYP3A5, PredDose, FirstDay, HCT) {
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
  Cp <- Cp * 1000
  Crbc = Cp * HCT * Bmax / (Kd + Cp) #concentration in red blood cells
  Cwb = (Cp + Crbc) #concentration in whole blood
  return(Cwb)
}

# Create OMEGA matric manually
omega <- matrix(data = c(0.14842,0.08379814,0.1380948,0,0,0,0.08379814,0.2558818,0,0,0,0,0.1380948,0,0.3342555,0,0,0,0,0,0,0.281337,0,0,0,0,0,0,0.05154826,0,0,0,0,0,0,0.891998), nrow = 6, ncol = 6)
colnames(omega) <- c("ECL", "EV1", "EQ", "EF", "EF_IOV", "EKa_IOV")
rownames(omega) <- c("ECL", "EV1", "EQ", "EF", "EF_IOV", "EKa_IOV")

# Create the tdmore model (based on the algebraic model)
tacro_algebraic <- tdmore(algebraic(tacrolimusAlgebraicFunction, output="Cwb"), omega = omega, res_var = list(errorModel("Cwb", prop=0.149)), iov=c("EKa_IOV", "EF_IOV"))

# Population plot
plot(tacro_algebraic, regimen, covariates=covariates, newdata=data.frame(TIME=seq(0, 7*24, by=1), Cwb=NA))

# Compute ipred based on observed data
ipred_algebraic <- estimate(tacro_algebraic, regimen=regimen, covariates=covariates, observed=observed)

# Find next dose
dose_algebraic <- findDose(ipred_algebraic, regimen=regimen, doseRows=c(4), target=data.frame(TIME=4*24, Cwb=13.5))
plot(ipred_algebraic, dose_algebraic$regimen, newdata=data.frame(TIME=seq(0, 7*24, by=1), Cwb=NA)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=12, ymax=15, fill="lightgreen", alpha=0.3)

options(scipen=8)
vcov(ipred_algebraic)

#______________________________________________________________________________________________
#----                              COMPARISON ODE VS ALGEBRAIC                             ----
#______________________________________________________________________________________________

expect_equal(tacro_ode$omega, tacro_algebraic$omega)
expect_equal(
  predict(tacro_ode, regimen=regimen, covariates=covariates,
        newdata=data.frame(TIME=seq(0, 50), Cwb=NA) )$Cwb,
  predict(tacro_algebraic, regimen=regimen, covariates=covariates,
        newdata=data.frame(TIME=seq(0, 50), Cwb=NA) )$Cwb
)
expect_equal(coef(ipred_ode), coef(ipred_algebraic), tolerance=5e-2)
expect_equal(vcov(ipred_ode), vcov(ipred_algebraic), tolerance=1e-1)
expect_equal(dose_algebraic$dose, dose_ode$dose, tolerance=1)
expect_equal(dose_algebraic$dose, 10.945, tolerance=1e-3)

#______________________________________________________________________________________________
#----                            TACROLIMUS ALGEBRAIC TIME-VARYING COVS                    ----
#______________________________________________________________________________________________
regimen <- data.frame(
  TIME=c(0,24,48),
  AMT=5,
  OCC=c(1,2,3)
)

# Example 1
covariatesVector=c(WT=70, HT=1.8, FEMALE=0, CYP3A5=0, PredDose=50, FirstDay=0, HCT=0.45)

covariates <- data.frame(TIME=c(0,24,48), EV1=c(1,0,2))
for(i in names(covariatesVector)) covariates[,i] <- covariatesVector[i]
example1 <- tacro_algebraic %>% predict(newdata=seq(0, 72, by = 0.1), regimen=regimen, covariates=covariates)
ggplot(example1, aes(x=TIME, y=Cwb)) + geom_line()

# Example 2
covariates <- data.frame(TIME=c(0,24,48), ECL=c(1,0,2))
for(i in names(covariatesVector)) covariates[,i] <- covariatesVector[i]
example2 <- tacro_algebraic %>% predict(newdata=seq(0, 72, by = 0.1), regimen=regimen, covariates=covariates)
ggplot(example2, aes(x=TIME, y=Cwb)) + geom_line()


