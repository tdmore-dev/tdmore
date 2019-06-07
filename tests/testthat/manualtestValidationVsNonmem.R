## Goal: to validate tdmore vs. Nonmem posthoc estimation
##
## Step 1: Generate a set of concentrations
## Step 2: Execute tdmore
## Step 3: execute nonmem on same model
## Step 4: compare results

# model from Delattre IK, Population pharmacokinetics of four β-lactams in critically ill septic patients comedicated with amikacin
# Meropenem model 1
library(RxODE)
library(tdmore)
library(dplyr)
library(ggplot2)

modelCode <- "
TVV1 = 24.4;
TVV2 = 7.01;
TVQ = 4.97;
TVCL = 9.87;
OV1 = 0.536;
OCL = 0.441;

EPS_PROP=0.371;

CL = TVCL * exp(ECL*OCL);
V1 = TVV1 * exp(EV1*OV1);
V2 = TVV2;
Q = TVQ;

K12 = Q/V1;
K21 = Q/V2;

d/dt(A1) = -CL/V1 * A1 - K12*A1 + K21 * A2;
d/dt(A2) = K12*A1 - K21 * A2;
CONC = A1 / V1;
"
rxModel <- RxODE::RxODE(modelCode)
mod <- rxModel %>% tdmore::tdmore(res_var=list(errorModel(prop=0.371)))
# 2g@0.5h as loading dose, 1g@0.5h at t=8
regimen <- data.frame(
  TIME=c(0, 8),
  AMT=c(2000, 1000),
  RATE=c(2000/0.5, 1000/0.5)
)
N <- 1000 #sample 100 patients
set.seed(0)
subjects <- data.frame(
  ID=1:N,
  ECL=rnorm(N),
  EV1=rnorm(N)
)
# Simulate from the model
db <- subjects %>% rowwise() %>% do({
  ETA <- .data
  myETA <- c(ECL=ETA$ECL, EV1=ETA$EV1)
  result <- mod %>%
    predict(
      newdata=data.frame(TIME=c(9, 16), CONC=NA),
      parameters=myETA,
      regimen=regimen)
  result$ID <- ETA$ID
  result
})
# Rich simulation profiles
dbrich <- subjects %>% rowwise() %>% do({
  ETA <- .data
  myETA <- c(ECL=ETA$ECL, EV1=ETA$EV1)
  result <- mod %>%
    predict(
      newdata=data.frame(TIME=seq(0, 16, by=0.5), CONC=NA),
      parameters=myETA,
      regimen=regimen)
  result$ID <- ETA$ID
  result
})
# Plot the simulated profiles, for error checking
ggplot(subset(dbrich, ID<=12)) + geom_line(aes(x=TIME, y=CONC)) + facet_wrap(~ID)
ggplot(dbrich %>% ddply(.(TIME), summarize, q50=median(CONC))) + geom_line(aes(x=TIME, y=q50))+
  geom_point(aes(x=TIME, y=CONC), data=db)

## For each row in DB, perform a TDM estimation
tdmorefits <- db %>% group_by(ID) %>% group_map( ~ {
  obs <- .x
  mod %>% estimate(observed=obs[,c("TIME", "CONC")], regimen=regimen)
})
estimationResults <- map_dfr(tdmorefits, function(fit) {
  fit %>% predict(newdata=data.frame(TIME=seq(0, 16, by=0.5), CONC=NA))
}, .id="ID")
coefResults <- map_dfr(tdmorefits, coef, .id="ID")

## Investigate a subject that was difficult to estimate
problematic <- mod %>% estimate(observed=db[db$ID==53,c("TIME", "CONC")], regimen=regimen, print.level=2)
problematicProfile <- problematic %>% profile(limits=c(-5, 5))
history <- data.frame(
  step=0:5,
  par1=c(0, 0.418, 0.894, 1.464, 2.300, 2.050),
  par2=c(0, -0.042, -0.289, -0.586, -1.024, -0.873),
  dpar1=-1*c(-4.180, -3.962, -4.583, -6.520, 11.816, -2.0845),
  dpar2=-1*c(0.417, 0.768, 1.242, 2.212, -4.307, 0.390)
)
## Should be used with 'logLik' data
scale_color <- function(x) {
  m <- mean(x)
  max <- max(x)
  min <- min(x)

  ifelse(x >= m,
                      # linear between median and upper limit
                      0.5 + 0.5*(exp(x)-exp(m)) / (exp(max) - exp(m)),
                      # loglinear between lower limit and median
                      0.5 * ( (x) - (m) ) / ( (min) - (m) )
                      )
}
problematicProfile$lik <- exp(problematicProfile$logLik)
ggplot(problematicProfile, aes(x=ECL, y=EV1)) +
  geom_tile(aes(fill=scale_color(logLik))) +
  #stat_contour(geom="polygon", aes(z=lik, fill=..level..), bins=50) +
  #stat_contour(geom="polygon", aes(z=exp(logLik), fill=..level..), bins=200) +
  geom_contour(aes(z=lik), bins=200) +
  geom_contour(aes(z=logLik), bins=5) +
  geom_point(data=data.frame(t(coef(problematic)))) +
  geom_segment(data=history, aes(x=par1, y=par2, xend=par1+dpar1, yend=par2+dpar2, group=step, color=factor(step)),
               arrow = arrow(length = unit(0.03, "npc")))+
  geom_text(data=history, aes(x=par1, y=par2, label=step, color=factor(step))) +
  scale_color_discrete() +
  scale_fill_gradient2(midpoint=0.5)

## Write a nonmem CSV file and control stream
nonmemdb <- ddply(db, .(ID), function(db) {
  nmregimen <- regimen
  nmregimen$EVID <- 1 # treatment
  nmregimen$ID <- db$ID[1]
  nmregimen$MDV <- 1

  db$EVID <- 0 #regular observation
  db$DV <- db$CONC
  db$MDV <- 0

  extraobs <- data.frame(
    ID=db$ID[1],
    EVID=2, #2 : do nothing, but give me the IPRED result at this time
    TIME=seq(0, 16, by=0.5),
    MDV=1
  )

  myDb <- rbind.fill(db, nmregimen, extraobs)
  myDb[order(myDb$ID, myDb$TIME, myDb$EVID),]
})
write.csv(nonmemdb, "nonmem.csv", row.names=F, na=".", quote=F)

modelCode <- "
$PROBLEM Posthoc estimation
$INPUT TIME CONC=DROP ID EVID DV MDV AMT RATE
$DATA nonmem.csv IGNORE=@
$SUBROUTINE ADVAN3 TRANS4
$PK
OV1 = 0.536;
OCL = 0.441;
CL = THETA(1) * EXP(ETA(1) * OCL)
V1 = THETA(2) * EXP(ETA(2) * OV1)
V2 = THETA(3)
Q = THETA(4)
S1=V1

$THETA
9.87 ;CL
24.4 ;V1
7.01 ;V2
4.97 ;Q

$OMEGA
;0.441 ; CL
;0.536 ; V1
1 ; CL
1 ; V1

$ERROR
IPRED=F
Y = F*(1+EPS(1))
;$SIGMA 0.371 ; Maybe we need the variance here as well??
$SIGMA 0.1376
;$EST METHOD=1 INTERACTION MAXEVAL=0 POSTHOC
$EST METHOD=1 INTERACTION MAXEVAL=0 POSTHOC LAPLACIAN
$TABLE ID TIME MDV DV EVID AMT RATE PRED IPRED ETA1 ETA2 FILE=sdtab ONEHEADER NOAPPEND
"
write(modelCode, "nonmem.ctl")
## execute nonmem code
psnOut <- shell("C:\\Perl64\\bin\\execute -clean=4 -nm_version=nm74g64 nonmem.ctl", intern=TRUE, mustWork=TRUE)
cat(psnOut)

sd <- read.table("sdtab", skip=1, header=TRUE)
ggplot(sd %>% subset(ID <= 12), aes(x=TIME)) +
  geom_point(aes(y=DV), data=function(db){subset(db, MDV==0)}) +
  geom_line(aes(y=IPRED)) +
  geom_line(aes(y=PRED), linetype=2) +
  facet_wrap(~ID)

## Compare the results
ggplot(sd %>% subset(ID <= 12), aes(x=TIME)) +
  geom_point(aes(y=DV), data=function(db){subset(db, MDV==0)}) +
  geom_line(aes(y=IPRED, color="nonmem"), size=1.5) +
  geom_line(aes(y=CONC, color="tdmore"), data=subset(estimationResults, ID<=12)) +
  facet_wrap(~ID)
comparison <- merge(coefResults, subset(sd, TIME==0.0 & EVID==1))
ggplot(comparison, aes(x=ECL, y=ETA1)) + geom_point() + geom_abline() + xlab("TDMore ETA_CL") + ylab("NONMEM ETA_CL") + geom_smooth(method="lm") + coord_fixed()
ggsave("compare ETA_CL.png", width=12, height=8)
ggplot(comparison, aes(x=EV1, y=ETA2)) + geom_point() + geom_abline() + xlab("TDMore ETA_V1") + ylab("NONMEM ETA_V1") + geom_smooth(method="lm") + coord_fixed()
ggsave("compare ETA_V1.png", width=12, height=8)

phi <- read.table("nonmem.phi", skip=1, header=TRUE)
vcovResults <- map_dfr(tdmorefits, function(x){
  result <- vcov(x)
  coef <- coef(x)
  c(coef, ETC.CL=result["ECL","ECL"], ETC.V1=result["EV1", "EV1"], ETC.CL.V1=result["ECL", "EV1"])
}, .id="ID")
vcovdb <- merge(phi, vcovResults, by.x="SUBJECT_NO", by.y="ID")

ggplot(vcovdb, aes(x=ETA.1., y=ECL)) + geom_point() + geom_abline() + labs(x="Nonmem estimate", y="TDMore estimate")
ggplot(vcovdb, aes(x=ETA.2., y=EV1)) + geom_point() + geom_abline() + labs(x="Nonmem estimate", y="TDMore estimate")
ggplot(vcovdb, aes(x=ETC.1.1., y=ETC.CL))+geom_point() + labs(x="Nonmem estimate", y="TDMore estimate") + geom_abline() + coord_fixed()
ggplot(vcovdb, aes(x=ETC.2.2., y=ETC.V1))+geom_point()
ggplot(vcovdb, aes(x=ETC.2.1., y=ETC.CL.V1))+geom_point()
