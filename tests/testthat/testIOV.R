library(tdmore)
library(nlmixr)
library(ggplot2)
library(testthat)
library(dplyr)

context("Test that IOV is working as expected")

set.seed(0)

regimen <- data.frame(
  TIME=c(0,24,48),
  AMT=150,
  OCC=c(1,2,2)
)

mod_1cpt_1 <- nlmixrUI(function(){
  ini({
    TVV <- 70
    TVKA <- 1
    TVCL <- 4
    ECL_IOV ~ 0.03
    ECL ~ 0.09
    EKA ~ 0.08
    EKA_IOV ~ 0.02
    EPS_PROP <- 0.1
  })
  model({
    CL <- TVCL * exp(ECL + ECL_IOV)
    V <- TVV
    KA <- TVKA * exp(EKA + EKA_IOV)
    d/dt(abs) = - abs*KA
    d/dt(center) = abs*KA - CL/V * center
    CONC = center / V
    CONC ~ prop(EPS_PROP)
  })
})

expect_error(mod_1cpt_1 %>% tdmore(iov="EKA_IOV2")) # Error is raised: 'IOV term(s) EKA_IOV2 not defined in model'
tdmore <- mod_1cpt_1 %>% tdmore(iov=c("EKA_IOV", "ECL_IOV"))

# Simple prediction test with IOV
data1 <- tdmore %>% predict(newdata=0:96, regimen=regimen, parameters=c(EKA_IOV=1, ECL_IOV=0, EKA_IOV=-2, ECL_IOV=0))
ggplot(data1, aes(x = TIME, y=CONC)) + geom_line()

observed <- data.frame(TIME=c(7, 23, 75), CONC=c(1.25, 0.55, 0.8))

# Simple estimate test with IOV
fit <- estimate(object = tdmore, regimen = regimen, observed = observed)
expectedValues <- c(ECL=-0.0102, EKA=0.0463, ECL_IOV=0.1026, EKA_IOV=0.0137, ECL_IOV=-0.1060, EKA_IOV=-0.0021)
expect_equal(round(coef(fit), digits=4), expectedValues)

data2 <- predict(fit, newdata=0:96, regimen=regimen)
ggplot(data2, aes(x = TIME, y=CONC)) + geom_line()

plot(fit, newdata=0:96)

# Check OMEGA matrix transformation
omega <- tdmore$omega
omega[2,1] <- 0.1
omega[1,2] <- 0.1
omega[3,1] <- 0.11
omega[1,3] <- 0.11
omega[4,1] <- 0.12
omega[1,4] <- 0.12
omega[3,2] <- 0.13
omega[2,3] <- 0.13
omega[4,2] <- 0.14
omega[2,4] <- 0.14
omega[4,3] <- 0.15
omega[3,4] <- 0.15

tdmore$omega <- omega
expandedOmega <- expandOmega(tdmore, 2)
expectedOmega <- matrix(c(0.09,0.13,0.1,0.14,0.1,0.14,0.13,0.08,0.11,0.15,0.11,0.15,0.1,0.11,0.03,0.12,0,0.12,0.14,0.15,0.12,0.02,0.12,0,0.1,0.11,0,0.12,0.03,0.12,0.14,0.15,0.12,0,0.12,0.02), nrow = 6, ncol = 6)
dimnames(expectedOmega) = list(names(expectedValues), names(expectedValues))
expect_equal(expandedOmega, expectedOmega)

# Test find dose
regimen <- data.frame(
  TIME=c(0,24,48,96),
  AMT=c(150,150,150,NA),
  OCC=c(1,2,2,3)
)

expect_error(findDose(fit, regimen, target=data.frame(TIME=120, CONC=1))) # Error: Number of occasions is different in tdmorefit regimen and findDose regimen

regimen <- data.frame(
  TIME=c(0,24,48,96),
  AMT=c(150,150,150,NA),
  OCC=c(1,2,2,2)
)

dose <- findDose(fit, regimen, target=data.frame(TIME=120, CONC=1), se.fit = T)
plot(fit, newdata=0:120, regimen = dose$regimen)
expect_equal(dose$dose[["dose.median"]], 206, tolerance=0.01)

# Test timevarying covariates
regimen <- data.frame(
  TIME=c(0,24,48),
  AMT=150,
  OCC=c(1,2,2)
)

mod_1cpt_1 <- nlmixrUI(function(){
  ini({
    TVV <- 70
    TVKA <- 1
    TVCL <- 4
    ECL_IOV ~ 0.03
    ECL ~ 0.09
    EKA ~ 0.08
    EKA_IOV ~ 0.02
    EPS_PROP <- 0.1
  })
  model({
    CL <- TVCL * (WT/70)^0.75 * exp(ECL + ECL_IOV)
    V <- TVV * WT/70
    KA <- TVKA * exp(EKA + EKA_IOV)
    d/dt(abs) = - abs*KA
    d/dt(center) = abs*KA - CL/V * center
    CONC = center / V
    CONC ~ prop(EPS_PROP)
  })
})
tdmore <- mod_1cpt_1 %>% tdmore(iov=c("EKA_IOV", "ECL_IOV"))
data3 <- tdmore %>% predict(newdata=seq(0, 96, by=0.1), regimen=regimen, parameters=c(EKA_IOV=1, ECL_IOV=0, EKA_IOV=-2, ECL_IOV=0), covariates=data.frame(TIME=c(0, 12.5, 24, 75), WT=c(70, 80, 85, 70)))

data4 <- rbind(data1 %>% mutate(type="No covariate"), data3 %>% mutate(type="WT covariate"))
ggplot(data4, aes(x = TIME, y=CONC, group=type, color=type)) + geom_line() +
  geom_line()


m1 <- tacrolimus_storset
regimen <- data.frame(
  TIME=seq(0, 6)*24,
  AMT=5,
  OCC=c(1,1,2,2,3,3,4)
)
covariates=c(WT=70, HT=1.8, FEMALE=0, CYP3A5=0, PredDose=50, FirstDay=0, HCT=0.45)
plot(m1, regimen, covariates=covariates, newdata=seq(0, 7*24, by=0.1) )

observed <- data.frame(TIME=3*24, Cwb=0.007)
ipred <- estimate(m1, regimen=regimen, covariates=covariates, observed=observed)
## Prediction at day 7 will be quite close to population, but prediction at day 4 should be more accurate
plot(ipred, newdata=seq(0, 7*24, by=0.1))
