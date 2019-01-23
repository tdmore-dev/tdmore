library(tdmore)
library(nlmixr)
library(ggplot2)
library(testthat)

context("Test that IOV is working as expected")

set.seed(0)

regimen <- data.frame(
  TIME=c(0,24,48),
  AMT=150,
  OCC=c(1,2,2)
)

regimen2 <- flatten(regimen)

mod_1cpt_1 <- nlmixrUI(function(){
  ini({
    TVV <- 70
    TVKA <- 1
    TVCL <- 4
    ECL ~ 0.09 # SD=0.3
    EKA ~ 0.09 # SD=0.3
    EKA_IOV ~ 0.09 # SD=0.3
    EPS_PROP <- 0.1
  })
  model({
    CL <- TVCL * exp(ECL)
    V <- TVV
    KA <- TVKA * exp(EKA + EKA_IOV)
    d/dt(abs) = - abs*KA
    d/dt(center) = abs*KA - CL/V * center
    CONC = center / V
    CONC ~ prop(EPS_PROP)
  })
})

expect_error(mod_1cpt_1 %>% tdmore(iov="EKA_IOV2")) # Error is raised: 'IOV term(s) EKA_IOV2 not defined in model'
m1 <- mod_1cpt_1 %>% tdmore(iov="EKA_IOV")


iov_parameters <- data.frame(EKA_IOV=c(0.1, 0.2))

debugonce(tdmore:::model_predict)
yep <- m1 %>% predict(newdata=0:96, regimen=regimen, iov_parameters=iov_parameters)

checkRegimen(regimen, "EKA_IOV")

getOccasionTimes(regimen)
