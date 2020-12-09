nlmixr::nlmixrUI(function(){
  ini({
    ## PK model
    TVCL <- 35
    TVVc <- 2500
    TVKa <- 0.125
    TVVp <- 750
    TVQ <- 1

    ECL ~ 0.06
    EVc ~ 0.06
    EKa ~ 2

    EPS_Prop <- 0.40

    ## PD model
    TVBASE <- 14
    TVKout <- 0.0003
    TVEMax <- 1
    TVEC50 <- 30
    TVKtol <- 0.00001

    EBASE ~ 0.8
    EKout ~ 0.5
    EEC50 ~ 3
    EKtol ~ 0.7

    EPS_Prop_SLD <- 0.15
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
