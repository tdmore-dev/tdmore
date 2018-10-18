library(tdmore)
library(RxODE)
library(nlmixr)
library(ggplot2)
library(testthat)

context("Test that the findDose method works as intended")

set.seed(0)

# Creating your model

modelCode <- function(){
  ini({
    TVV1 <- 24.4;
    TVV2 <- 7.01;
    TVQ <- 4.97;
    TVCL <- 9.87;
    ECL ~ 0.194 # This value corresponds to OMEGA_CL (44% SD)
    EV1 ~ 0.287 # This value corresponds to OMEGA_V1 (54% SD)
    EPS_PROP <- 0.371 # Proportional error (37% SD)
  })
  model({
    CL <- TVCL * exp(ECL)
    V1 <- TVV1 * exp(EV1)
    V2 <- TVV2
    Q <- TVQ
    K12 <- Q/V1
    K21 <- Q/V2

    d/dt(center) = - CL/V1 * center - K12*center + K21 * periph
    d/dt(periph) = K12*center - K21 * periph

    CONC = center / V1
    CONC ~ prop(EPS_PROP) # Proportional error linked to the PK model
  })
}

nlmixrUI <- nlmixrUI(modelCode)
tdmore <- tdmore(nlmixrUI)

# Predicting new data

regimen <- data.frame(
  TIME=c(0, 8, 16),            # Every 8 hour and for 1 day, an injection is given
  AMT=c(1000, 1000, 1000),     # 1g is administered
  RATE=c(1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
)

data <- predict(
  object = tdmore,
  newdata = data.frame(TIME = seq(0, 24, by = 0.5), CONC = NA),
  regimen = regimen,
  se = TRUE
)


ggplot(data, aes(x=TIME, y=CONC)) +
  geom_ribbon(aes(fill="Population", ymin=CONC.lower, ymax=CONC.upper), fill="steelblue2", alpha=0.15) +
  geom_line(aes(color="Population"), data=data) +
  scale_color_manual(values=c("steelblue2"))

# Estimating individual parameters

pred <- estimate(tdmore = tdmore, regimen = regimen)

observed <- data.frame(TIME=c(9, 16), CONC=c(30, 8))

ipred <- estimate(tdmore = tdmore, observed = observed, regimen = regimen)

data <- predict(ipred, newdata=data.frame(TIME=seq(0, 24, 0.1), CONC=NA), se=TRUE)

ggplot(data, aes(x=TIME))  +
  geom_line(aes(color="Individual", y=CONC.median)) +
  geom_ribbon(aes(fill="Individual", ymin=CONC.lower, ymax=CONC.upper), fill="tomato1", alpha=0.10) +
  geom_line(aes(color="Population", y=CONC), data=predict(pred, newdata=seq(0, 24, 0.1))) +
  geom_point(aes(y=CONC), data=observed) +
  scale_color_manual(values=c("tomato1", "steelblue2"))

# Finding the right dose to give

newRegimen <- data.frame(
  TIME=c(0, 8, 16, 24),              # A fourth dose on the second day is added
  AMT=c(1000, 1000, 1000, NA),       # Adding an unknown dose on the second day
  RATE=c(1000, 1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
)

recommendation1 <- findDose(
  ipred,
  regimen = newRegimen,
  interval = c(100, 5000),
  target = data.frame(TIME = 32, CONC = 3.1)
)

expect_equal(round(recommendation1$dose, digits=2), 296.62)

# Finding the right dose to give (with 95% CI)

recommendation2 <- findDose(
  ipred,
  regimen = newRegimen,
  interval = c(100, 5000),
  target = data.frame(TIME = 32, CONC = 3.1),
  se.fit = T,
  level = 0.95,
  mc.maxpts = 50,
  extendInt = "yes"
)
expect_equal(round(recommendation2$dose$dose.median, digits=2), 224.83)

# Continue
regimen <- recommendation2$regimen
observed <- data.frame(TIME=c(9, 16, 20, 32), CONC=c(30, 8, 15, 3.1+0.2))

ipred <- estimate(tdmore = tdmore, observed = observed, regimen = regimen)

plot(ipred)
