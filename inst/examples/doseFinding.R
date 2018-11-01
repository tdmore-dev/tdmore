## Lots of different examples will be needed:
# simple dose optimization to reach a next target asap
# to reach target in steady-state
# with loading dose
# with PD target
# with target window
rm(list=ls(all=TRUE))
library(RxODE)
library(ggplot2)
library(tdmore)

# Population pharmacokinetics of tacrolimus in
# adult kidney transplant recipients, Staatz et al 2002
modelCode <- "
CL = 23.6 * exp(ETA1*0.42);
Vc = 1070 * exp(ETA2*1.11);
ka=4.48;
CONC = centr / Vc * 1000;
d/dt(abs) = -ka*abs;
d/dt(centr) = ka*abs - CL/Vc*centr;
"
model <- RxODE::RxODE(modelCode) %>%
  tdmore(res_var=list(errorModel(var="CONC", add=3.7))) #Model has 3.7 ng/mL additive error

regimen <- data.frame(
  TIME=c(0, 24),
  AMT=c(15, 15),
  II=c(0, 24),
  ADDL=c(0, 10)
)
observed <- data.frame(TIME=c(2.4, 23), CONC=c(10, 5))
ipredfit <- model %>%
  estimate(observed, regimen)

D <- findDose(ipredfit, regimen=regimen, target=data.frame(TIME=48, CONC=13.5))

z <- plot(ipredfit)


gridExtra::grid.table(observed)

pred <- model %>% estimate(regimen=regimen)
predmc <- predict(pred, newdata = data.frame(TIME=seq(0, 120, length.out=500), CONC=NA), se.fit=TRUE, level=0.95, mc.maxpts=500, .progress="text")
ggplot(predmc, aes(x=TIME, color="Population prediction", fill="Population prediction")) +
  geom_ribbon(aes(ymin=CONC.lower, ymax=CONC.upper), alpha=0.3) +
  geom_line(aes(y=CONC.median)) +
  xlab("Time (hour)") +
  ylab("Concentration (ng/mL)") +
  coord_cartesian(ylim=c(0, 50)) +
  scale_x_continuous(breaks=seq(0, 120, by=24))

ipredfit <- model %>%
  estimate(observed, regimen)

melt <- function(x, se=FALSE) {
  measure.vars <- colnames(x)
  measure.vars <- measure.vars[measure.vars != "TIME"]
  vars <- measure.vars
  if( se ) {
    for(i in c(".upper", ".lower")) vars <-c(vars, paste0(measure.vars, i))
  }
  tmp <- reshape::melt(x, id.vars="TIME")
  if(se) {
    result <- subset(tmp, variable %in% measure.vars)
    result$value.upper <- subset(tmp, variable %in% paste0(measure.vars, ".upper"))$value
    result$value.lower <- subset(tmp, variable %in% paste0(measure.vars, ".lower"))$value
  } else {
    result <- tmp
  }
  result
}

newdata = data.frame(TIME=seq(0, 100, length.out=500), CONC=NA)
tdmorefit <- ipredfit




ipred <- tdmorefit %>% predict(newdata, regimen=regimen) %>% melt
ggplot(ipred, aes(x=TIME, y=value)) + geom_line() + geom_point(aes(y=CONC), data=observed)

ipredre <- tdmorefit %>% predict(newdata, regimen=regimen, se.fit=TRUE, level=0.95, mc.maxpts = 500, .progress="text") %>% melt(se=TRUE)
pred <- estimate(tdmorefit$tdmore, regimen=regimen) %>% predict(newdata) %>% melt
predre <- estimate(tdmorefit$tdmore, regimen=regimen) %>% predict(newdata, se.fit=TRUE, level=0.95, mc.maxpts=500, .progress="text") %>% melt(se=TRUE)
obs <- model.frame(tdmorefit, se=TRUE, level=0.95) %>% melt(se=TRUE)

ggplot(mapping=aes(x=TIME, y=value)) +
  geom_ribbon(aes(fill="Population (95% CI)", ymin=value.lower, ymax=value.upper), data=predre, alpha=0.3)+
  geom_ribbon(aes(fill="Fit (95% CI)", ymin=value.lower, ymax=value.upper), data=ipredre, alpha=0.3)+
  geom_line(aes(color="Population"), data=pred) +
  geom_line(aes(color="Fit"), data=ipred) +
  geom_point(data=observed, aes(x=TIME, y=CONC)) +
  geom_errorbar(aes(ymax=value.upper, ymin=value.lower), data=obs) +
  xlab("Time (hour)") + ylab("Tac concentration (ng/mL)") + coord_cartesian(ylim=c(0, 50)) +
  scale_x_continuous(breaks=seq(0, 120, by=24))








tdmorefit <- ipredfit
newregimen <- data.frame(
  TIME=c(0, 24),
  AMT=c(15, 10),
  II=c(0, 24),
  ADDL=c(0, 10)
)
ipred <- tdmorefit %>% predict(newdata, regimen=newregimen) %>% melt
ipredre <- tdmorefit %>% predict(newdata, regimen=newregimen, se.fit=TRUE, level=0.95, mc.maxpts = 500, .progress="text") %>% melt(se=TRUE)
pred <- estimate(tdmorefit$tdmore, regimen=newregimen) %>% predict(newdata) %>% melt
predre <- estimate(tdmorefit$tdmore, regimen=newregimen) %>% predict(newdata, se.fit=TRUE, level=0.95, mc.maxpts=500, .progress="text") %>% melt(se=TRUE)
obs <- model.frame(tdmorefit, se=TRUE, level=0.95) %>% melt(se=TRUE)

ggplot(mapping=aes(x=TIME, y=value)) +
  geom_ribbon(ymin=12, ymax=15, aes(fill="Therapeutic window"), alpha=0.2, data=data.frame(TIME=c(0, 100), value=10)) +
  geom_ribbon(aes(fill="Population (95% CI)", ymin=value.lower, ymax=value.upper), data=predre, alpha=0.3)+
  geom_ribbon(aes(fill="Fit (95% CI)", ymin=value.lower, ymax=value.upper), data=ipredre, alpha=0.3)+
  geom_line(aes(color="Population"), data=pred) +
  geom_line(aes(color="Fit"), data=ipred) +
  geom_point(data=observed, aes(x=TIME, y=CONC)) +
  #geom_errorbar(aes(ymax=value.upper, ymin=value.lower), data=obs) +
  xlab("Time (hour)") + ylab("Tac concentration (ng/mL)") + coord_cartesian(ylim=c(0, 30))




z <- plot(ipredfit, newdata=data.frame(TIME=seq(0, 100, length.out=200), CONC=NA), se.obs=FALSE, maxpts=500, .progress="text")
print(z)

searchspace <- searchspace() %>%
  combinations.searchspace(AMT=c(1, 2, 5, 10, 15))

evaluators <- list()
evaluators[["C0_next > 0.05"]] = function(tdmorefit, regimen) {
  tmax <- max(tdmorefit$observed$TIME)
  myReg <- flatten(regimen)
  nextRegimen <- match(TRUE, myReg$TIME > tmax)
  nextRegimenTrough <- myReg$TIME[ nextRegimen+1]
  pred <- predict(tdmorefit, newdata=nextRegimenTrough-0.01, regimen)
  pred$CONC > 0.05
}

evaluators[["DosePrice"]] = function(tdmorefit, regimen) {
  totalDose <- sum( regimen$AMT * (1 + regimen$ADDL) )
  totalDose * 150 #150 euro / mg
}

# x <- evaluate(ipred, regimen, searchspace, evaluators)
# plot.possibilities(ipred, regimen, searchspace, evaluators)
