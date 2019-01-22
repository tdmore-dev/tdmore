library(tdmore)
library(nlmixr)
library(ggplot2)
library(testthat)

context("Test that IOV is working as expected")

set.seed(0)

regimen <- data.frame(
  TIME=c(0,24),
  AMT=150
)

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
