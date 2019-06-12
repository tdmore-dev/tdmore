##
## Script name: validation.R
##
## Purpose of script:
## Numerical validation of tdmore vs NONMEM
##
## Author: Ruben Faelens
##
## Date Created: Wed Nov 07 12:02:54 2018
##
## Copyright (c) Ruben Faelens, 2018
## Email: ruben.faelens@gmail.com
##
## ---------------------------
##
## Notes:
##     A huge thank you to D. Kang for sharing the nonmem code used in the paper below:
##     Kang, Dongwoo, et al. "Standard error of empirical bayes estimate in NONMEM® VI." The Korean Journal of Physiology & Pharmacology 16.2 (2012): 97-106.
##
##     This was the basis for confirming NONMEM predictions match the R theoretical predictions.
## ---------------------------
nmPath <- "C:/nm74g64/util/"


library(magrittr)

wd <- tempfile(pattern="tdmore-validation")
dir.create(wd) %>% stopifnot()



# NONMEM estimation results -----------------------------------------------
## Copy all NONMEM data
theoPath <- file.path(nmPath, "THEO")
file.copy(theoPath, wd) %>% stopifnot()

phenoPath <- file.path(nmPath, "PHENO")
file.copy(phenoPath, wd) %>% stopifnot()

## Copy all NONMEM control streams
nmCtlFiles <- list.files(pattern="*.CTL", path="inst/validation/NONMEM", full.names=T)
for(f in nmCtlFiles)
  file.copy(f, wd) %>% stopifnot()

## Execute all NONMEM control streams
nmfePattern <- if(.Platform$OS.type=="windows") "^nmfe.*\\.bat$" else "^nmfe.*"
nmfe <- list.files(pattern=nmfePattern, path=nmPath, full.names=T)[1]


originalWd <- getwd()
setwd(wd)
nmCtlFiles <- list.files(pattern="^.*\\.CTL$", path=wd, full.names=T)
results <- plyr::alply(nmCtlFiles, 1, function(f){
  command <- paste(nmfe, f, paste0(f, ".lst"))
  system(command, intern=TRUE)
}, .progress="text")
setwd(originalWd)

## Read all the results
eta <- list.files(pattern="^.*\\.phi", path=wd, full.names=T)
db <- plyr::adply(eta, 1, function(f){
  db <- read.table(f, header=T, skip=1, check.names=F)
  db$file <- basename(f)
  db
})
db$type <- factor(
  substr(db$file, 0, 1),
  levels=c("B", "C"),
  labels=c("THEO", "PHENO")
)
db$method <- substr(db$file, 4, nchar(db$file)-4)
db <- db[, !colnames(db) %in% c("X1", "SUBJECT_NO", "file")]
results_nonmem <- db

## Plot the results
ggplot(results_nonmem, aes(x=ID, y=`ETA(1)`, color=method)) + geom_point() + facet_wrap(~type, scales="free")



# Helper functions to convert tdmore to nonmem.phi format -----------------
asPhi <- function(estimate) {
  par <- coef(estimate)
  df1 <- as.data.frame( t(par) )
  colnames(df1) <- paste0("ETA(", seq_along(par), ")")

  vcov <- vcov(estimate)
  n <- nrow(vcov)
  result <- c()
  for(i in 1:n) {
    for(j in 1:i) {
      value <- vcov[i, j]
      names(value) <- paste0("ETC(",i,",",j,")")
      result <- c(result, value)
    }
  }
  df2 <- data.frame( t(result), check.names = F)

  cbind(df1, df2)
}

# TDMore estimation results: Theophylline -----------------------------------------------
## From B11FOCEi.ext
## THETA1       THETA2       THETA3       SIGMA(1,1)   SIGMA(2,1)   SIGMA(2,2)   OMEGA(1,1)   OMEGA(2,1)   OMEGA(2,2)   OMEGA(3,1)   OMEGA(3,2)   OMEGA(3,3)
## 1.49072E+00  3.24470E+01  8.72901E-02  1.70386E-02  0.00000E+00  8.28440E-02  4.34717E-01  5.69202E-02  1.93995E-02 -6.48548E-03  1.20830E-02  2.02334E-02
## TODO: Estimate the model parameters using nlmixr?
## TODO: Can we increase atol and rtol for RxODE?

modelCode <- function() {
  ini({
    #THETA1 <- 2
    #THETA2 <- 50
    #THETA3 <- 0.1
    THETA1 <- 1.49072E+00
    THETA2 <- 3.24470E+01
    THETA3 <- 8.72901E-02
    #variances
    #                     OMEGA(1,1)   OMEGA(2,1)   OMEGA(2,2)   OMEGA(3,1)   OMEGA(3,2)   OMEGA(3,3)
    ETA1 + ETA2 + ETA3 ~ c(4.34717E-01,  5.69202E-02,  1.93995E-02, -6.48548E-03,  1.20830E-02,  2.02334E-02)
    EPS1 <- 0.1305320 #1.70386E-02 variance
    EPS2 <- 0.2878263 #8.28440E-02 variance
  })
  model({
    KA=THETA1*exp(ETA1)
    V=THETA2*exp(ETA2)
    K=THETA3*exp(ETA3)

    d/dt(abs) = -KA*abs
    d/dt(centr) = KA*abs - K*centr
    #TODO: we do not support algebraic equations?
    #CONC=DOSE / V * KA/(KA-K)*(exp(-K*t) - exp(-KA*t))
    CONC=centr/V
    CONC ~ prop(EPS1) + add(EPS2)
  })
}

m1 <- nlmixr::nlmixrUI(modelCode)
m1tdm <- tdmore(m1, atol=1e-12, rtol=1e-12) ## We need high precision to match the estimates from the algebraic implementation in nonmem

regimen <- data.frame(TIME=0, AMT=320) ## fixed

input <- read.table(file=theoPath, na=".")
names(input) <- c("ID", "AMT", "TIME", "DV", "BWT")

results_m1 <- input %>% group_by(ID) %>% do({
  observed <- .data
  obs <- observed[, c("TIME", "DV")]
  colnames(obs) <- c("TIME", "CONC")
  est <- estimate(m1tdm, obs, regimen)

  # tibble(ipred=list(est))
  cbind(asPhi(est),
       data.frame(
         ID=observed$ID[1],
         OBJ=est$ofv,
         type="THEO",
         method="tdmore"
      ))
})

# TDMore estimation results: Theophylline -----------------------------------------------
## From C11FOCEi.ext
#ITERATION    THETA1       THETA2       THETA3       SIGMA(1,1)   SIGMA(2,1)   SIGMA(2,2)   OMEGA(1,1)   OMEGA(2,1)   OMEGA(2,2)   OBJ
#33           4.47972E-01  1.38052E+00  4.83244E-03  6.32087E-03  0.00000E+00  3.89899E+00  1.68867E-02  8.72258E-03  1.41515E-02    578.27818681292047
## TODO: Estimate the model parameters using nlmixr?
## TODO: Can we increase atol and rtol for RxODE?

modelCode <- function() {
  ini({
    THETA1 <- 4.47972E-01
    THETA2 <- 1.38052E+00
    THETA3 <- 4.83244E-03
    #variances
    ETA1 + ETA2 ~ c(1.68867E-02,  8.72258E-03,  1.41515E-02 )
    EPS1 <- 0.0795039
    EPS2 <- 1.9745860
  })
  model({
    TVV  = THETA1 + (BWT/1.5)**THETA2
    TVK  = THETA3
    V    = TVV * exp(ETA1)
    K    = TVK * exp(ETA2)

    d/dt(centr) = -K*centr
    #TODO: we do not support algebraic equations?
    CONC=centr/V
    CONC ~ prop(EPS1) + add(EPS2)
  })
}

m2 <- nlmixrUI(modelCode)
m2tdm <- tdmore(m2)

input <- read.table(file=phenoPath, na=".", nrows=744) ## PHENO has 'FIN' as the last line...
names(input) <- c("ID", "TIME", "AMT", "BWT", "APGR", "DV", "MDV", "EVID")

results_m2 <- input %>% group_by(ID) %>% do({
  input <- .data
  observed <- filter(input, EVID==0 & MDV==0)
  obs <- observed[, c("TIME", "DV")]
  colnames(obs) <- c("TIME", "CONC")

  regimen <- filter(input, EVID==1)
  regimen <- regimen[, c("TIME", "AMT")]

  covariates <- input[, c("TIME", "BWT")]

  est <- estimate(m2tdm, obs, regimen, covariates)

  cbind(asPhi(est),
        data.frame(
          ID=observed$ID[1],
          OBJ=est$ofv,
          type="PHENO",
          method="tdmore"
        ))
})


# Compare results ---------------------------------------------------------

fulldb <- plyr::rbind.fill(results_nonmem, results_m1) %>% plyr::rbind.fill(results_m2) %>% filter(method %in% c("tdmore", "FOCEi")) %>%
  tidyr::gather(key=variable, value=value, -ID, -method, -type, -OBJ) %>% dplyr::rename(meth=method)
ggplot(fulldb, aes(x=ID, y=value, color=meth)) + geom_line(na.rm=TRUE) + facet_grid(variable~type, scales="free")
ggsave("eta_profiles.png", width=12, height=8)

comparison <- fulldb %>%
  dplyr::select(-OBJ) %>%
  tidyr::spread(meth, value) %>% ## cast has a bug when using the 'method' column...
  dplyr::mutate(delta=tdmore-FOCEi)
ymax <- abs(max(comparison$delta, na.rm=T)) * 5

ggplot(comparison, aes(x=ID, y=delta, size=abs(delta))) +
  geom_point(na.rm=TRUE) +
  facet_grid(type~variable, scales="free_y") +
  scale_y_continuous(labels=scales::percent, limits=c(-ymax, ymax)) +
  scale_x_continuous(breaks=seq(1, 100, by=2)) +
  geom_hline(yintercept=0) +
  coord_flip() +
  labs(title="TDMore vs FOCEi", x="delta", size="delta")
ggsave("eta_profiles_delta.png", width=16, height=9)
