library(tdmore)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
set.seed(1234)

# Model taken from literature: Soulele, K., et al.
# 'Population pharmacokinetics of fluticasone propionate/salmeterol using two different dry powder inhalers.'
# European Journal of Pharmaceutical Sciences 80 (2015): 33-42."
myModel <- nlmixr::nlmixrUI(function(){
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

m1 <- tdmore(myModel)

regimen <- rbind(
  data.frame(TIME=0, AMT=500), #loading dose
  data.frame( #maintenance dose
    TIME=sort( c(8+(0:20)*24, 20+(0:20)*24) ), #at 08:00 and 20:00
    AMT=88
  )
)
obs <- data.frame(
  TIME=(3*24+12)+(0:2)*7*24, #observe every 3rd day of the week at noon
  CONC=NA
)

# predict generates N samples of a tdmorefit with uncertainty
# If we use the population fit, then this will represent a sample across the population.
# I.e. N virtual patients
# In a model with covariates, you should either sample the covariates from an external database,
# or generate them from a distribution.
N <- 10
population <- estimate(m1, regimen=regimen)
dbPred <- predict(population, newdata=obs, se=T, level=NA, mc.maxpts=N) %>%
  rename(ID=sample) %>%
  select(ID, TIME, CONC)
# We add residual error by sampling from the original model
# This is probably too pessimistic; residual error is rarely truly random...
db <- dbPred %>%
  model.frame(m1, data=., se=TRUE, level=NA) #sample residual error

ggplot(dbPred, aes(x=TIME, y=CONC)) +
  geom_point(alpha=0.3) +
  geom_line(aes(group=ID), data=db) +
  geom_smooth()

posthocFit <- posthoc(m1, regimen=regimen, observed=db, se.fit=F, control=list(factr=1e14))
posthocFit %>%
  mutate(coef = map(fit, ~tibble::enframe(coef(.x)))) %>%
  tidyr::unnest(coef) %>%
  ggplot(aes(x=ID, y=value)) +
    geom_point() +
    facet_wrap(~name) +
    labs(title="Eta estimates for posthoc fit")

ggplot(db, aes(x=TIME, y=CONC)) +
  geom_point() +
  geom_line(data= posthocFit %>% mutate(ipred = map(fit, predict)) %>% tidyr::unnest(ipred) ) +
  facet_wrap(~ID) +
  labs(title="Individual fits using all available data")

## Input data: 'm1', the tdmore model
## Input data: 'db' with columns ID, TIME and CONC
## Input data: 'regimen' with columns TIME and AMT.
## Note that with true retrospective data, the regimen would be different per individual
prosevalDb <- proseval(m1, regimen=regimen, observed=db, control=list(factr=1e14))

# Calculate predictions for all data, using the fits
prosevalPred <- prosevalDb %>%
  left_join(db %>% tidyr::nest(-ID), by="ID") %>%
  mutate(ipred = pmap(list(fit, newdata=data), predict) )

target <- prosevalPred %>%
  tidyr::unnest(data, ipred, .sep=".") %>%
  group_by(ID, OBS) %>%
  filter(row_number() == OBS[1]+1) %>%  #pick the next week as a target
  mutate(IRES = data.CONC-ipred.CONC)

ggplot(target, aes(x=IRES)) +
  stat_ecdf() +
  geom_vline(xintercept=0) +
  facet_wrap(~OBS) +
  labs(title="Absolute prediction error on day i+1\nUsing on all data available on day i")

parameters <- prosevalDb %>%
  group_by(ID, OBS) %>%
  group_modify(~predict(.x$fit[[1]], newdata=0)) %>%
  tidyr::gather(key=key, value=value, Ka:V2) %>%
  group_by(ID, key) %>%
  mutate(relative=(value-value[1])/value)
parameters %>%
  ggplot(aes(x=OBS, y=relative*100)) +
  facet_wrap(~key, scales="free") +
  geom_point() +
  geom_line(aes(group=ID)) +
  labs(title="Evolution of parameter values over time",
       y="Relative change (%)")

ggplot(target, aes(x=ipred.CONC / data.CONC)) +
  stat_ecdf() +
  scale_x_log10(limits=c(0.5, 2)) +
  geom_vline(xintercept=1) +
  geom_vline(xintercept=c(0.81, 1.22), linetype=2) +
  geom_label(aes(x=1.8, y=0.5, label=label), data=. %>% group_by(OBS) %>% summarize(label=paste0( round(mean(between(ipred.CONC/data.CONC, 0.81, 1.22))*100), "%"))) +
  facet_wrap(~OBS) +
  labs(title="Relative prediction error on day i+1\nUsing on all data available on day i")
