library(tdmore)
nlmixr::nlmixrUI(function(){
  ini({
    TVKa <- 3
    TVCL <- 700    # L/h
    TVV1 <- 60000  # L
    TVV2 <- 6000   # L
    TVQ <- 250     # L/h

    EKa ~ 0.04
    ECL ~ 0.15
    EV1 ~ 0.09
    EV2 ~ 0.20
    EQ ~ 0.10

    EPS_ADD <- 2
    EPS_PROP <- 0.12
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
}) %>%
  tdmore::tdmore() %>%
  tdmore::metadata(tdmore::output(name="CONC", label="Concentration", unit="mg/L")) %>%
  tdmore::metadata(tdmore::formulation(name="Drug", unit="mg", dosing_interval=8, default_value=5, round_function=function(x){round(x*2)/2})) %>%
  tdmore::metadata(tdmore::target(min=12, max=20))
