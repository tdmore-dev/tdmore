# This project
`tdmore` provides an easy interface to execute post-hoc bayesian estimation of individual profiles, and to find the best dose to put the patient on target.

The package is intended to make it easy to define your own models, and easily create a dose decision support tool for physicians.

# License
This project does not include a license. This means that all work is under exclusive copyright. Nobody else can use, copy, distribute or modify this work.

The Github terms of service apply. We allow others to view and fork the repository. Please note that this is not sufficient to then copy, distribute or modify this work further.

Please see https://choosealicense.com/no-permission/ for more information.

Through publishing, we allow others to use this R package and to perform dose adaptation. Installing this package using `devtools::install_github` is allowed.

The official copyright holder of this work is the KU Leuven university.

# Limitations
This software is a research project, and cannot be considered as a medical device. It is not a substitute for clinical reasoning.

# How to install
```R
devtools::install_github("tdmore-dev/tdmore")
```

We suggest you also install `tidyverse` and `RxODE`. The latter requires a working C and fortran compiler to work.

# How to use
This example code uses `RxODE` to define the model. 

```R
library(RxODE)
library(ggplot2)
library(tdmore)

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
