# Khosravan, Reza, et al. "Population pharmacokinetic/pharmacodynamic modeling of sunitinib by dosing schedule in patients with advanced renal cell carcinoma or gastrointestinal stromal tumor." Clinical pharmacokinetics 55.10 (2016): 1251-1269.
nlmixr::nlmixrUI(function(){
  ini({
    ## PK model sunitinib
    TVCL <- 34.1
    TVVc <- 2700
    TVKa <- 0.126
    TVVp <- 774
    TVQ <- 0.688

    ECL ~ 0.060516 # 24.6%
    EVc ~ 0.052900 # 23.0%
    EKa ~ 2.755600 # 166%

    EPS_Prop <- 0.417

    ## PD model SLD (Tumor sum of longest diameters)
    TVBASE <- 14.3 #cm
    TVKout <- 0.000267
    TVEMax <- 1
    TVEC50 <- 30.5
    TVKtol <- 0.0000141

    EBASE ~ 0.840889 # 91.7%
    EKout ~ 0.521284 # 72.2%
    EEC50 ~ 3.312400 # 182%
    EKtol ~ 0.720801 # 84.9%

    EPS_Prop_SLD <- 0.143
  })
  model({
    CL <- TVCL * exp(ECL)
    Vc <- TVVc * exp(EVc)
    Vp <- TVVp
    Q <- TVQ
    K12 <- Q/Vc
    K21 <- Q/Vp
    Ke <- CL/Vc
    Ka <- TVKa*exp(EKa)
    BASE <- TVBASE * exp(EBASE)
    Kout <- TVKout * exp(EKout)
    EMax <- TVEMax
    EC50 <- TVEC50 * exp(EEC50)
    Ktol <- TVKtol * exp(EKtol)
    Kin <- Kout*BASE

    SLD(0) = BASE;

    d/dt(depot) = -Ka*depot
    d/dt(center) = Ka*depot - Ke*center - K12*center + K21*periph
    d/dt(periph) = K12*center - K21*periph
    DRUG = EMax*(center/Vc) / (EC50 + center/Vc)
    d/dt(SLD) = Kin*(1-DRUG) - Kout*SLD*(1+exp(-Ktol*t))

    CONC = center/Vc
    CONC ~ prop(EPS_Prop) | center

    SLD ~ prop(EPS_Prop_SLD) | SLD
  })
})
