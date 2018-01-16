# This project
`tdmore` provides an easy interface to execute post-hoc bayesian estimation of individual profiles, and to find the best dose to put the patient on target.

The package is intended to make it easy to define your own models, and easily create a dose decision support tool for physicians.

# How to install
```R
devtools::install_github("rfaelens/tdmore")
```

We suggest you also install `tidyverse` and `RxODE`. The latter requires a working C and fortran compiler to work.

# How to use
This example code uses `RxODE` to define the model. 

```R
library(RxODE)
library(ggplot2)

modelCode <- "
CL = 3.7 * exp(ETA1*0.19);
Vc = 61 * exp(ETA2*0.28);
ka=3.7;
CONC = centr / Vc;
d/dt(abs) = -ka*abs;
d/dt(centr) = ka*abs - CL/Vc*centr;
"
model <- RxODE::RxODE(modelCode) %>%
  tdmore(prop=0.23) #Model has 23% proportional error

regimen <- data.frame(
  TIME=seq(0, 7)*24,
  AMT=5 #5mg
)
observed <- data.frame(TIME=2, CONC=0.04)

ipred <- model %>%
  estimate(observed, regimen)

print(summary(ipred))

z <- ggplot(ipred %>% profile(maxpts=20), aes(x=ETA1, y=ETA2, z=logLik)) + geom_contour()
print(z)

newdata = data.frame(TIME=seq(0, 12, length.out=50), CONC=NA)
z <- plot.tdmorefit(ipred, newdata, .progress="text")
print(z)
```
