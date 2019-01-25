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
    ECL_IOV ~ 0.02
    ECL ~ 0.09 # SD=0.3
    EKA ~ 0.09 # SD=0.3
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


#debugonce(tdmore:::model_predict)

data1 <- tdmore %>% predict(newdata=0:96, regimen=regimen, parameters=c(EKA_IOV=1, ECL_IOV=0, EKA_IOV=-2, ECL_IOV=0))
ggplot(data1, aes(x = TIME, y=CONC)) + geom_line()

#observed <- data.frame(TIME=c(20, 75), CONC=c(0.8, 0.9))
observed <- data.frame(TIME=c(7, 23, 75), CONC=c(1.25, 0.55, 0.8))

#debugonce(tdmore:::estimate)
#debug(tdmore:::pop_ll)
fit <- estimate(object = tdmore, regimen = regimen, observed = observed)

#debugonce(tdmore:::predict.tdmorefit)
data2 <- predict(fit, newdata=0:96, regimen=regimen)
ggplot(data2, aes(x = TIME, y=CONC)) + geom_line()

#debugonce(tdmore:::plot.tdmorefit)
plot(fit, newdata=0:96)
