# Vignettes  {#vignettes}

## Running TDMore

### What's TDMore? {-}

TDMore is an R package that attemps to provide a global framework for making Therapeutic Drug Monitoring software. It is designed not only for the pharmacometricians but also for the physicians.

### Creating your model {-}

There are several ways telling TDMore what your model is. The easiest way is to provide the model in the form of a nxlmir model. The following example shows how a 2-compartment-PK model can be written. For the purpose of this demonstration, the 'Meropenem' PK model will be used. As Meropenem is given by injection into a vein, the absorption compartment is omitted.



```r
library(nlmixr)

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
```

Once the model is created, a TDMore object can be instantiated as follows:


```r
library(tdmore)

nlmixrUI <- nlmixrUI(modelCode)
tdmore <- tdmore(nlmixrUI)
```

The TDMore object is now ready for use.

### Predicting new data {-}

TDMore can be used to run simulations, based on the model defined in the previous step. For doing so, the regimen first needs to be specified. In the case of Meropenem, an 30-min injection is given into the central compartment every 8 hour. This can be written as follows:


```r
regimen <- data.frame(
  TIME=c(0, 8, 16),            # Every 8 hour and for 1 day, an injection is given
  AMT=c(1000, 1000, 1000),     # 1g is administered
  RATE=c(1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
)
```

Let's now simulate the population PK model for 1 day. This can be done using the TDMore predict() function. A dataframe with all the times to predict (and respective NA concentration) is given to the 'newdata' argument. Field 'se' is set to true, it tells TDMore we are interested to simulate the model with residual variability.


```r
data <- predict(
  object = tdmore,
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
```

<img src="Vignettes_files/figure-html/population_prediction-1.png" width="768" style="display: block; margin: auto;" />

The plot above shows the median of the population (typical profile) and its 95% confidence interval. The range of the confidence interval can be changed using the argument 'level'.

### Estimating individual parameters {-}

This section will show you how the individual parameters can be estimated, based on some observed data. Let's first estimate the parameters of a typical individual. This is achieved by calling the estimate() function.


```r
pred <- estimate(tdmore = tdmore, regimen = regimen)
coef(pred)
```

```
## ECL EV1 
##   0   0
```

Both eta's ECL and EV1 have been estimated to 0. This is not surprising, as zero eta's best describe the population average. We can also look at the uncertainty on these estimates.


```r
vcov(pred)
```

```
##       ECL   EV1
## ECL 0.194 0.000
## EV1 0.000 0.287
```

This uncertainty is equal to the population inter-individual variability (OMEGA matrix). Now, let's assume blood samples have been collected for a subject X at different times. For example, blood samples were collected at times 9h and 16h on the first day. This can be translated in TDMore as follows (note that the concentrations are purely fictive):


```r
observed <- data.frame(TIME=c(9, 15), CONC=c(30, 2))
```

We can ask TDMore to re-estimate the parameters for this specific individual:


```r
ipred <- estimate(tdmore = tdmore, observed = observed, regimen = regimen)
coef(ipred)
```

```
##        ECL        EV1 
##  0.2597598 -0.1320256
```

Eta's obtained by calling 'coef' on pred maximise altogether the likelihood for this specific subject. The variance-covariance matrix shows the uncertainty of the individual estimates, and their correlation.


```r
vcov(ipred)
```

```
##            ECL        EV1
## ECL 0.03154084 0.04680264
## EV1 0.04680264 0.12759969
```

Predictions for the population (pred) and this specific subject (ipred) can be compared using the following snippet: 


```r
data <- predict(ipred, newdata=data.frame(TIME=seq(0, 24, 0.1), CONC=NA), se=TRUE)
ggplot(data, aes(x=TIME))  +
  geom_line(aes(color="Individual", y=CONC.median)) +
  geom_ribbon(aes(fill="Individual", ymin=CONC.lower, ymax=CONC.upper), fill="tomato1", alpha=0.10) +
  geom_line(aes(color="Population", y=CONC), data=predict(pred, newdata=seq(0, 24, 0.1))) +
  geom_point(aes(y=CONC), data=observed) +
  scale_color_manual(values=c("tomato1", "steelblue2")) +
  scale_y_log10()
```

<img src="Vignettes_files/figure-html/pred_ipred_predictions-1.png" width="768" style="display: block; margin: auto;" />

Ipred and pred are shown respectively in red and blue. A 95% confidence interval has been added around the ipred prediction. Note that the default TDMore plotting function can also be used to obtain the exact same plot in a single line of code:


```r
plot(ipred, newdata=data.frame(TIME=seq(0.1, 24, by=0.1), CONC=NA)) + scale_y_log10()
```

### Finding the right dose to give {-}

A very interesting feature in TDMore is the possibility to ask the framework the next dose to be given knowing all the previous observations that were collected and some known end-points. For example, assume we still collected the same two observations on the first day, we would like to find the best first dose to be given on the second day. We would like to reach the trough concentration of 3.10 mg/L as much as possible. This can be expressed, as follows:


```r
newRegimen <- data.frame(
  TIME=c(0, 8, 16, 24),              # A fourth dose on the second day is added
  AMT=c(1000, 1000, 1000, NA),       # Adding an unknown dose on the second day
  RATE=c(1000, 1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
)

recommendation <- findDose(
  ipred,
  regimen = newRegimen,
  interval = c(100, 5000),
  target = data.frame(TIME = 32, CONC = 8)
  )
summary(recommendation)
```

```
## $dose
## [1] 4540.072
## 
## $regimen
##   TIME      AMT RATE
## 1    0 1000.000 2000
## 2    8 1000.000 2000
## 3   16 1000.000 2000
## 4   24 4540.072 2000
```

The result of the findDose() routine is shown above. It tells us that XXX mg (approximately) is the recommended starting dose on the second day. The following code helps up verify this visually.


```r
# Recommended regimen for individual can be directly accessed from the recommendation object
ipredRecommendedRegimen <- recommendation$regimen

# Population regimen can be updated using the 'updateRegimen' method, a 4th dose of 1000 is used
predUpdatedRegimen <- updateRegimen(regimen = newRegimen, newDose = 4500)
print(predUpdatedRegimen) # Check pred regimen
```

```
##   TIME  AMT RATE
## 1    0 1000 2000
## 2    8 1000 2000
## 3   16 1000 2000
## 4   24 4500 2000
```


```r
ipred <- estimate(tdmore = tdmore, observed = observed, regimen = ipredRecommendedRegimen)
pred <- estimate(tdmore = tdmore, regimen = predUpdatedRegimen)
```


```r
plot(ipred, newdata=seq(0, 32, by=0.1)) + geom_hline(yintercept=8) + scale_y_log10()
```

<img src="Vignettes_files/figure-html/find_dose_pred_ipred_predictions-1.png" width="768" style="display: block; margin: auto;" />

The plot above demonstrates that the individual is reaching the trough concentration quite well after the first administration on the second day. 

## Multiple endpoints

### Writing and testing the PK model {-}

In this vignette, we will learn how TDMore can deal with multiple endpoints. The PK/PD models chosen to illustrate this section are based on the following paper: "Population PK/PD modeling of Sunitinib by dosing schedule in patients with advanced renal cell carcinoma or gastrointestinal stromal tumor.". PDF file is available here. Let's start writing the Sunitinib PK model.  
  



```r
library(nlmixr)

modelCode <- function(){
  ini({
    TVCL <- 34.1
    TVVc <- 2700
    TVKa <- 0.126
    TVVp <- 774
    TVQ <- 0.688

    ECL ~ 0.060516 # 24.6%
    EVc ~ 0.052900 # 23.0%
    EKa ~ 2.755600 # 166%

    EPS_Prop <- 0.417
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

    d/dt(depot) = -Ka*depot
    d/dt(center) = Ka*depot - Ke*center - K12*center + K21*periph
    d/dt(periph) = K12*center - K21*periph

    CONC = center/Vc
    CONC ~ prop(EPS_Prop)
  })
}
```

The TDMore object is instantiated as follows:


```r
library(tdmore)

nlmixrUI <- nlmixrUI(modelCode)
tdmore <- tdmore(nlmixrUI)
```

A basic regimen can be created to test that the model is properly running. 50 mg Sunitinib is given for a week.


```r
regimen <- data.frame(
  TIME=0,   # First dose time: t=0h
  AMT=50,   # Dose amount: 50 mg
  II=24,    # Dose interval: 24h
  ADDL=7  # Additional doses: 7
)

times <- seq(0, 1*7*24, by=1) # Observation times
```

This regimen can be plotted using the default TDMore plotting function. It shows the typical value of the population and the between-subject variability (95% confidence interval).


```r
plot(tdmore, regimen, newdata=times)
```

<img src="Vignettes_files/figure-html/sunitinib_pk_model-1.png" width="768" style="display: block; margin: auto;" />

### Adding a PD model {-}

Many different PD models are described in the paper mentionned above. We will focus on the simplest one: a PD model related to the target tumor's sum of the longest diameters (abbreviated SLD), which corresponds to the main efficacy endpoint of Sunitinib. Adding this PD model to the existing PK model is done as follows. Please note the mandatory '|' nlmixr syntax used to describe the residual variability of both endpoints (CONC and SLD).


```r
modelCode <- function(){
  ini({
    # PK model sunitinib
    TVCL <- 34.1
    TVVc <- 2700
    TVKa <- 0.126
    TVVp <- 774
    TVQ <- 0.688

    ECL ~ 0.060516 # 24.6%
    EVc ~ 0.052900 # 23.0%
    EKa ~ 2.755600 # 166%

    EPS_Prop <- 0.417 # Proportional error 1 (related to CONC)

    # PD model SLD (Tumor sum of longest diameters)
    TVBASE <- 14.3 #cm
    TVKout <- 0.000267
    TVEMax <- 1
    TVEC50 <- 30.5
    TVKtol <- 0.0000141

    EBASE ~ 0.840889 # 91.7%
    EKout ~ 0.521284 # 72.2%
    EEC50 ~ 3.312400 # 182%
    EKtol ~ 0.720801 # 84.9%

    EPS_Prop_SLD <- 0.143 # Proportional error 2 (related to SLD)
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
    d/dt(SLD) = Kin*(1-DRUG) - Kout*SLD*(1+exp(-Ktol*t)) # PD model ode

    CONC = center/Vc
    CONC ~ prop(EPS_Prop) | center # Define error model 1

    SLD(0) = BASE                  # Set SLD initial value at time t=0
    SLD ~ prop(EPS_Prop_SLD) | SLD # Define error model 2
  })
}
nlmixrUI <- nlmixrUI(modelCode)
tdmore <- tdmore(nlmixrUI, maxsteps=1E3*500) # Old tdmore object is overridden
```

Let's now have a look at the evolution of SLD over time. To have a good overview, we'll observe SLD for 40 weeks. The default plot shows once again the typical value and the between-subject variability (95% CI).


```r
regimen <- data.frame(
  TIME=0,
  AMT=50,
  II=24,
  ADDL=40*7
)
times <- seq(0, 40*7*24, by=1)

# Note that 'SLD' is specified in the newdata dataframe to select the right endpoint
# If not specified, both endpoints are printed
plot(tdmore, regimen, newdata=data.frame(TIME=times, SLD=NA))
```

<img src="Vignettes_files/figure-html/sunitinib_pd_model-1.png" width="768" style="display: block; margin: auto;" />

### Estimating individual parameters {-}

Assume SLD is measured at week 0 and week 30 for a certain individual. Model parameters can be estimated and visualised as follows:

```r
observed <- data.frame(TIME=c(0, 30*7*24), CONC=NA, SLD=c(25, 14))
ipred <- estimate(tdmore = tdmore, observed = observed, regimen = regimen)
plot(ipred, newdata=data.frame(TIME=times, SLD=NA))
```

<img src="Vignettes_files/figure-html/sunitinib_pd_ipred_sld-1.png" width="768" style="display: block; margin: auto;" />

Coefficients show that the PD model could be estimated correctly. 

```r
coef(ipred)
```

```
##           ECL           EVc           EKa         EBASE         EKout 
##  6.762903e-05  9.882637e-07 -7.529206e-07  5.396754e-01 -6.433377e-02 
##         EEC50         EKtol 
##  3.814871e-03  1.165208e-02
```

```r
vcov(ipred)
```

```
##                 ECL           EVc           EKa         EBASE
## ECL    6.051170e-02 -3.202960e-08  7.977322e-08 -3.281933e-05
## EVc   -3.202960e-08  5.290001e-02 -1.574797e-08 -1.274422e-07
## EKa    7.977322e-08 -1.574797e-08  2.755604e+00  5.372490e-07
## EBASE -3.281933e-05 -1.274422e-07  5.372490e-07  1.253910e-02
## EKout  1.947017e-04  5.265617e-07 -1.905303e-06  3.441526e-02
## EEC50 -2.383770e-04 -8.312489e-07  3.327081e-06 -1.808524e-03
## EKtol -2.901414e-05 -1.094317e-07  4.564505e-07 -5.593560e-03
##               EKout         EEC50         EKtol
## ECL    1.947017e-04 -2.383770e-04 -2.901414e-05
## EVc    5.265617e-07 -8.312489e-07 -1.094317e-07
## EKa   -1.905303e-06  3.327081e-06  4.564505e-07
## EBASE  3.441526e-02 -1.808524e-03 -5.593560e-03
## EKout  3.849601e-01  1.071095e-02  3.278036e-02
## EEC50  1.071095e-02  3.299272e+00 -1.598567e-03
## EKtol  3.278036e-02 -1.598567e-03  7.239781e-01
```

However, TDMore is not able to estimated PK parameters because no data was provided, as confirmed by the following plot (IPRED strictly equal to PRED, only the last week is shown).


```r
plot(ipred, newdata=data.frame(TIME=seq(29*7*24, 30*7*24, by=1), CONC=NA), se.fit=F)
```

<img src="Vignettes_files/figure-html/sunitinib_pd_ipred_conc-1.png" width="768" style="display: block; margin: auto;" />
