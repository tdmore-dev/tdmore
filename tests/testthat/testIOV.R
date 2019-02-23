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

expect_error(mod_1cpt_1 %>% tdmore(iov="EKA_IOV2"),
             regexp = "IOV term.*EKA_IOV2 not defined in model") # Error is raised: 'IOV term(s) EKA_IOV2 not defined in model'
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
#omega[2,1] <- 0.1
#omega[1,2] <- 0.1
#omega[3,1] <- 0.11
#omega[1,3] <- 0.11
omega[4,1] <- 0.12
omega[1,4] <- 0.12
omega[3,2] <- 0.13
omega[2,3] <- 0.13
#omega[4,2] <- 0.14
#omega[2,4] <- 0.14
#omega[4,3] <- 0.15
#omega[3,4] <- 0.15

tdmore$omega <- omega
expandedOmega <- expandOmega(tdmore, 2)
expectedOmega <- matrix(c(0.09, 0.13, 0, 0, 0, 0, 0.13, 0.08, 0, 0, 0, 0, 0,
                          0, 0.03, 0.12, 0, 0, 0, 0, 0.12, 0.02, 0, 0, 0, 0, 0, 0, 0.03,
                          0.12, 0, 0, 0, 0, 0.12, 0.02), nrow = 6, ncol = 6)
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
    ECL ~ 0.09
    EKA ~ 0.08
    ECL_IOV ~ 0.03
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

data4Check <- data4 %>% filter(TIME==55) %>% select(CONC)
expect_equal(data4Check, data.frame(CONC=c(1.717376, 1.460061)), tolerance=1e-6)

tdmore$omega[1,2] <- tdmore$omega[2,1] <- sqrt(0.09)*sqrt(0.08)*-0.1 #IIV correlation
tdmore$omega[3,4] <- tdmore$omega[4,3] <- sqrt(0.03)*sqrt(0.02)*0.3 #IOV correlation
tdmore$omega

test_that("expandOmega works as intended", {
  iovIndex <- rownames(tdmore$omega) %in% tdmore$iov
  expect_equal(
    expandOmega(tdmore, 0),
    tdmore$omega[!iovIndex, !iovIndex]
  )
  expect_equal(
    expandOmega(tdmore, 1),
    tdmore$omega
  )

  expectedResult <- Matrix::bdiag(
    tdmore$omega[!iovIndex, !iovIndex],
    tdmore$omega[iovIndex, iovIndex],
    tdmore$omega[iovIndex, iovIndex]
  ) %>% as.matrix()
  myNames <- rownames(tdmore$omega)
  myNames <- c(myNames[!iovIndex], rep(myNames[iovIndex], 2) )
  rownames( expectedResult ) <- myNames
  colnames( expectedResult ) <- myNames
  expect_equal(
    expandOmega(tdmore, 2),
    expectedResult
  )
})
