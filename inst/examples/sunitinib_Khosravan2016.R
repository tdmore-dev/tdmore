library(nlmixr)
## Model based on literature
## Khosravan, Reza, et al. "Population pharmacokinetic/pharmacodynamic modeling of sunitinib by dosing schedule in patients with advanced renal cell carcinoma or gastrointestinal stromal tumor." Clinical pharmacokinetics 55.10 (2016): 1251-1269.
modelCode <- function(){
  ini({
    ## PK model sunitinib
    TVCL=34.1
    TVVc=2700
    TVKa=0.126
    TVTLag=0.527
    TVVp=774
    TVQ=0.688
    TVAGE_CL=-0.00702
    TVRAC_CL=-0.152
    TVSEX_CL=-0.193
    TVTUMR_CL=0.293
    TVBWT_Vc=0.281
    TVSEX_Vc=-0.213
    TVTUMR_Vc=0.420

    EPS_Prop = 0.417

    ECL ~ 0.060516 #24.6%
    EVc ~ 0.052900 #23.0%
    EKa ~ 2.755600 #166%

    ## PD model SLD
    TVBASE = 14.3 #cm
    TVKout=0.000267
    TVEMax = 1
    TVEC50=30.5
    TVKtol=0.0000141
    TVBEC_BASE=0.574
    TVRAC_BASE=-0.348
    TVSCH_BASE=-0.430
    TVSCH_Kout=1.01
    TVSCH_EC50=2.43
    TVTUMR_EC50=4.82

    EBASE ~ 0.840889# 91.7%
    EKout ~ 0.521284#72.2%
    EEC50 ~ 3.312400#182%
    EKtol ~ 0.720801#84.9%

    EPS_Prop_SLD = 0.143
  })
  model({
    CL <- TVCL * exp(ECL)
    Vc <- TVVc * exp(EVc)
    Vp <- TVVp
    Q <- TVQ
    K12 <- Q/Vc
    K21 <- Q/Vp
    K0 <- CL/Vc
    Ka <- TVKa*exp(EKa)

    BASE=TVBASE * exp(EBASE)
    KOut = TVKout*exp(EKout)
    EMax = TVEMax
    EC50=TVEC50*exp(EEC50)
    Ktol = TVKtol*exp(EKtol)
    CONC = center / V1
    DRUG = EMax * CONC / (EC50 + CONC)
    Kin = Kout*BASE

    SLD(0) = BASE


    d/dt(abs) = -Ka*abs
    d/dt(center) = Ka*abs - K0 * center - K12*center + K21 * periph
    d/dt(periph) = K12*center - K21 * periph
    d/dt(SLD) = Kin*(1-DRUG) - Kout*SLD*exp(-Ktol*t)

    ## nlmixr does something weird for the residual error model
    ## It seems you cannot concurrently specify two different RE outputs?
    CONC ~ prop(EPS_PROP) | DVTYPE==1
    SLD ~ prop(EPS_Prop_SLD) | DVTYPE==2
  })
}
nlmixrModel <- nlmixrUI(modelCode)

library(tdmore)
m1 <- tdmore(nlmixrModel)

regimen <- data.frame(
  TIME=0,
  AMT=50,
  II=24,
  ADDL=4*7
)

data <- predict(
object = m1,
newdata = data.frame(TIME = seq(0, 24, by = 0.5), CONC = NA),
regimen = regimen,
se = TRUE
)

library(ggplot2)
ggplot(data, aes(x=TIME, y=CONC)) +
geom_ribbon(aes(fill="Population", ymin=CONC.lower, ymax=CONC.upper), fill="steelblue2", alpha=0.15) +
geom_line(aes(color="Population"), data=data) +
scale_color_manual(values=c("steelblue2")) +
scale_y_log10()

observed <- data.frame(TIME=c(9, 15), CONC=c(30, 2), SLD=c(27,25))

ipred <- estimate(tdmore = m1, observed = observed, regimen = regimen)
coef(ipred)
vcov(ipred)

data <- predict(ipred, newdata=data.frame(TIME=seq(0, 24, 0.1), CONC=NA), se=TRUE)
ggplot(data, aes(x=TIME))  +
  geom_line(aes(color="Individual", y=CONC.median)) +
  geom_ribbon(aes(fill="Individual", ymin=CONC.lower, ymax=CONC.upper), fill="tomato1", alpha=0.10) +
  geom_line(aes(color="Population", y=CONC), data=predict(pred, newdata=seq(0, 24, 0.1))) +
  geom_point(aes(y=CONC), data=observed) +
  scale_color_manual(values=c("tomato1", "steelblue2")) +
  scale_y_log10()
