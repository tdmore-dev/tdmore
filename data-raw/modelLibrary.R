
default_model <- nlmixrUI(function(){
  ini({
    TVKA <- 3.7
    TVV1 <- 61
    TVQ <- 10
    TVCL <- 3.7
    ECL ~ 0.0784 #ETA1, 28%
    EV1 ~ 0.0361 #ETA2, 19%
    EPS_PROP <- 0.23 # Proportional error, 23% SD
  })
  model({
    KA <- TVKA
    CL <- TVCL * (WT/70)^0.75 * exp(ECL)
    V1 <- TVV1 * exp(EV1)
    V2 <- V1
    Q <- TVQ
    K12 <- Q/V1
    K21 <- Q/V2

    d/dt(depot) = -KA*depot
    d/dt(center) = KA*depot - CL/V1 * center - K12*center + K21 * periph
    d/dt(periph) = K12*center - K21 * periph

    CONC = center / V1
    CONC ~ prop(EPS_PROP)
  })
})
usethis::use_data(default_model, overwrite=TRUE)

two_error_models <- nlmixrUI(function(){
  ini({
    TVKA <- 3.7
    TVV1 <- 61
    TVQ <- 10
    TVCL <- 3.7
    ECL ~ 0.0784 #ETA1, 28%
    EV1 ~ 0.0361 #ETA2, 19%
    EPS_PROP <- 0.001 # Proportional error, 10% SD
    EPS_PROP2 <- 0.001 # Proportional error, 10% SD
  })
  model({
    KA <- TVKA
    CL <- TVCL * (WT/70)^0.75 * exp(ECL)
    V1 <- TVV1 * exp(EV1)
    V2 <- V1
    Q <- TVQ
    K12 <- Q/V1
    K21 <- Q/V2

    d/dt(depot) = -KA*depot
    d/dt(center) = KA*depot - CL/V1 * center - K12*center + K21 * periph
    d/dt(periph) = K12*center - K21 * periph

    CONC_C = center / V1
    CONC_C ~ prop(EPS_PROP) | center

    CONC_P = periph / V2
    CONC_P ~ prop(EPS_PROP2) | periph
  })
})
usethis::use_data(two_error_models, overwrite=TRUE)

meropenem_model <- nlmixrUI(function(){
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
})
usethis::use_data(meropenem_model, overwrite=TRUE)

meropenem_model_wt <- nlmixrUI(function(){
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
    CL <- TVCL * (WT/70)^0.75 * exp(ECL)
    V1 <- TVV1 * (WT/70) * exp(EV1)
    V2 <- TVV2 * (WT/70)
    Q <- TVQ * (WT/70)^0.75
    K12 <- Q/V1
    K21 <- Q/V2

    d/dt(center) = - CL/V1 * center - K12*center + K21 * periph
    d/dt(periph) = K12*center - K21 * periph

    CONC = center / V1
    CONC ~ prop(EPS_PROP) # Proportional error linked to the PK model
  })
})
usethis::use_data(meropenem_model_wt, overwrite=TRUE)

meropenem_omega0_model <- nlmixrUI(function(){
  ini({
    TVV1 <- 24.4;
    TVV2 <- 7.01;
    TVQ <- 4.97;
    TVCL <- 9.87;
    ECL ~ 0.0 # This value corresponds to OMEGA_CL (0% SD)
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
})
usethis::use_data(meropenem_omega0_model, overwrite=TRUE)

meropenem_1param_model <- nlmixrUI(function(){
  ini({
    TVV1 <- 24.4;
    TVV2 <- 7.01;
    TVQ <- 4.97;
    TVCL <- 9.87;
    ECL <- 0.0 # Fixed to 0. Not considered as a parameter in TDMore.
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
})
usethis::use_data(meropenem_1param_model, overwrite=TRUE)

# Model taken from literature: Soulele, K., et al.
# 'Population pharmacokinetics of fluticasone propionate/salmeterol using two different dry powder inhalers.'
# European Journal of Pharmaceutical Sciences 80 (2015): 33-42."

fluticasone_model <- nlmixrUI(function(){
  ini({
    TVKa <- 3.87
    TVCL <- 659    # L/h
    TVV1 <- 56900  # L
    TVV2 <- 5550   # L
    TVQ <- 259     # L/h

    EKa ~ 0.04507129  # 0.2123**2
    ECL ~ 0.1535856   # 0.3919**2
    EV1 ~ 0.09223369  # 0.3037**2
    EV2 ~ 0.208301    # 0.4564**2
    EQ ~ 0.1015697    # 0.3187**2

    EPS_ADD <- 1.91
    EPS_PROP <- 0.117
  })
  model({
    Ka <- TVKa * exp(EKa)
    CL <- TVCL * exp(ECL)
    V1 <- TVV1 * exp(EV1)
    V2 <- TVV2 * exp(EV2)
    Q <- TVQ * exp(EQ)
    K12 <- Q/V1
    K21 <- Q/V2

    d/dt(center) = - CL/V1 * center - K12*center + K21 * periph
    d/dt(periph) = K12*center - K21 * periph

    CONC = center / V1 * 1000
    CONC ~ prop(EPS_PROP) + add(EPS_ADD)
  })
})
usethis::use_data(fluticasone_model, overwrite=TRUE)

# Khosravan, Reza, et al. "Population pharmacokinetic/pharmacodynamic modeling of sunitinib by dosing schedule in patients with advanced renal cell carcinoma or gastrointestinal stromal tumor." Clinical pharmacokinetics 55.10 (2016): 1251-1269.
sunitinib_pkpd_model <- nlmixrUI(function(){
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

    d/dt(depot) = -Ka*depot
    d/dt(center) = Ka*depot - Ke*center - K12*center + K21*periph
    d/dt(periph) = K12*center - K21*periph
    DRUG = EMax*(center/Vc) / (EC50 + center/Vc)
    d/dt(SLD) = Kin*(1-DRUG) - Kout*SLD*(1+exp(-Ktol*t))

    CONC = center/Vc
    CONC ~ prop(EPS_Prop) | center

    SLD(0) = BASE
    SLD ~ prop(EPS_Prop_SLD) | SLD
  })
})
usethis::use_data(sunitinib_pkpd_model, overwrite=TRUE)
