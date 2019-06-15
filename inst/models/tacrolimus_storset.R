# Størset, Elisabet, et al. "Improved prediction of tacrolimus concentrations early after kidney transplantation using theory‐based pharmacokinetic modelling." British journal of clinical pharmacology 78.3 (2014): 509-523.
# Fat free mass was determined by:
#Janmahasatian, Sarayut, et al. "Quantification of lean bodyweight." Clinical pharmacokinetics 44.10 (2005): 1051-1065.
#36/205 for CYP3A=0, CYP3A=1
tdmore::tdmore( nlmixr::nlmixrUI(function(){
  ini({
    TVKa <- 1.01
    #TVTLag <- 0.41  #to be included in model definition itself...

    TVCLp = 811
    TVCLp_CYP3A5 = 1.30
    TVV1p = 6290
    TVQp=1200
    TVV2p=32100

    TVF_PredEMax = 0.67 #relative
    TVF_Pred50 = 35 #mg
    TVF_FirstDay = 2.68 #relative
    TVF_CYP3A5  = 0.82

    # var <- list(CL=0.14842, V1=0.2558818, Q=0.3342555)
    #40%, 54%, 63% resp; 0.43 corrCL/V1, 0.62 corrCL/Q
    ECL + EV1 + EQ ~ c(0.14842,
                       0.08379814, 0.2558818,
                       0.1380948,  0,         0.3342555)
    # Corresponding variance calculated using CV = sqrt( exp(omega^2) -1 )
    # omega^2 = log( CV^2 + 1 )
    EF ~ 0.281337 #57%

    EF_IOV ~ 0.05154826 #23%
    EKa_IOV ~ 0.891998 #120%

    EPS_PROP = 0.149 #standard deviation
  })
  model({
    BMI = WT / HT^2 ;
    FFM = 9.27E3 * WT / (6.68E3 + 216*BMI);
    if(FEMALE) FFM=9.27E3 * WT / (8.78E3 + 244*BMI);
    Ka <- TVKa*exp(EKa_IOV)

    CLp <- TVCLp * (FFM/60)^0.75 * TVCLp_CYP3A5^CYP3A5 * exp(ECL)
    V1p <- TVV1p * exp(EV1)
    V2p <- TVV2p
    Qp <- TVQp * exp(EQ)

    K12 <- Qp/V1p
    K21 <- Qp/V2p
    Ke <- CLp/V1p

    #Covariates: PredDose, FirstDay, CYP3A5
    TVF = 1 * (1 - TVF_PredEMax*PredDose / (PredDose + TVF_Pred50)) * TVF_FirstDay^FirstDay * TVF_CYP3A5^CYP3A5
    F = TVF * exp(EF + EF_IOV)

    d/dt(abs) = -Ka*abs
    d/dt(center) = F*Ka*abs - Ke*center - K12*center + K21*periph
    d/dt(periph) = K12*center - K21*periph

    Bmax = 418 #ug/L erythrocytes
    Kd=3.8 #ug/L plasma
    Cp = center/V1p * 1000 #free concentration in plasma, in ng/mL or ug/L
    Crbc = Cp * HCT * Bmax / (Kd + Cp)
    Cwb = (Cp + Crbc) #in ng/mL
    Cwb ~ prop(EPS_PROP)
  })
}), iov=c("EF_IOV", "EKa_IOV"))
