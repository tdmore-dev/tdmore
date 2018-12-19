# THEOPP from nm74
theopp <- read.table("data-raw/THEOPP", na='.')
colnames(theopp) <- c("ID", "AMT", "TIME", "CONC", "WT")
theopp$EVID <- ifelse(is.na(theopp$AMT), 0, 1)
#theopp$WT <- zoo::na.locf(theopp$WT)  #not strictly required?
usethis::use_data(theopp, overwrite = TRUE)

# Estimate using nlmixr
library(nlmixr)
library(tidyverse)
modelCode <- function() {
  ini({
    LTHETA1 <- log(2)
    LTHETA2 <- log(50)
    LTHETA3 <- log(0.1)
    #variances
    ETA1 ~ 1
    ETA2 ~ 2
    ETA3 ~ 1
    #ETA1 + ETA2 + ETA3 ~ c(0.2, 0.0, 0.2, 0.0, 0.0, 0.2)
    EPS1 <- 0.1
    EPS2 <- 0.1
  })
  model({
    KA=exp(LTHETA1+ETA1)
    V=exp(LTHETA2+ETA2)
    K=exp(LTHETA3+ETA3)

    d/dt(abs) = -KA*abs
    d/dt(centr) = KA*abs - K*centr
    CONC=centr/V
    CONC ~ prop(EPS1) + add(EPS2)
  })
}


result <- theopp %>% rename(DV=CONC) %>% nmDataConvert %>%
  nlmixr(modelCode, data=., est="saem")
theopp_nlmixr <- result
usethis::use_data(theopp_nlmixr, overwrite = TRUE)
