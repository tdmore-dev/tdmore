
# PHENO from nm74
phenoRaw <- readLines("data-raw/PHENO")
n <- which(grepl(phenoRaw, pattern="FIN$")) #file ends with a line "FIN"
pheno <- read.table("data-raw/PHENO", nrows = n-1, na='.')
colnames(pheno) <- c("ID", "TIME", "AMT", "WT", "APGR", "DV", "EVID", "MDV")
usethis::use_data(pheno, overwrite = TRUE)



# Estimate using nlmixr
library(nlmixr)
library(tidyverse)
modelCode <- function() {
  ini({
    LTHETA2 <- log(50)
    LTHETA3 <- log(0.1)
    #variances
    ETA2 ~ 2
    ETA3 ~ 1
    #ETA1 + ETA2 + ETA3 ~ c(0.2, 0.0, 0.2, 0.0, 0.0, 0.2)
    EPS2 <- 0.1
  })
  model({
    LWT = log(WT / 70)
    V=exp(LTHETA2+LWT+ETA2)
    CL=exp(LTHETA3+LWT*0.75+ETA3)

    d/dt(centr) = -CL/V*centr
    CONC=centr/V
    CONC ~ add(EPS2)
  })
}


result <- nlmixr(modelCode, data=pheno, est="focei")
pheno_nlmixr <- result
usethis::use_data(pheno_nlmixr, overwrite = TRUE)
