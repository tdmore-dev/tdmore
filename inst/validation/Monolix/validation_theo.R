##
## Script name: validation.R
##
## Purpose of script:
## Numerical validation of tdmore vs Monolix
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
## This script verifies whether
## a) The structural model results in the same predicted values
## b) The mode of the monolix individual estimates matches tdmore
## c) Whether IOV and time-varying covariates handling works the same way
##
## ---------------------------

library(testthat)


# Load the Monolix Theophylline demo --------------------------------------
library(tidyverse)
library(MlxConnectors)

MlxConnectors::initializeMlxConnectors("monolix")
mlxInfo <- MlxConnectors::getMlxEnvInfo()[[1]]
mlxDirectory <- get('MLX_DIRECTORY', attr(mlxInfo, 'envir') )

demoPath <- file.path(mlxDirectory, "resources/demos/monolix/1.creating_and_using_models/1.1.libraries_of_models")

MlxConnectors::loadProject(file.path(demoPath, "theophylline_project.mlxtran") )
MlxConnectors::getContinuousObservationModel()
MlxConnectors::getObservationInformation()
#MlxConnectors::setErrorModel(CONC="combined2")
data <- MlxConnectors::getData()

MlxConnectors::runPopulationParameterEstimation()
MlxConnectors::runConditionalModeEstimation()

MlxConnectors::getEstimatedRandomEffects()

mlxEta <- MlxConnectors::getEstimatedRandomEffects(method="conditionalMode")
mlxEta <- mlxEta$conditionalMode %>% mutate(id=as.numeric(as.character(id)))
obs <- MlxConnectors::getObservationInformation()$CONC

mlxPar <- MlxConnectors::getEstimatedIndividualParameters(method="conditionalMode")
mlxPred <- MlxConnectors::computePredictions(individualParameters = mlxPar$conditionalMode)$Cc
mlxPred <- cbind( obs, data.frame(IPRED=mlxPred))


# Run tdmore, using the parameter estimates from Monolix ------------------
library(tdmore)
m1 <- algebraic(fun=function(t, TIME, AMT, ka_pop, V_pop, Cl_pop, eta_ka, eta_V, eta_Cl) {
  ka = ka_pop * exp(eta_ka)
  V = V_pop * exp(eta_V)
  Cl = Cl_pop * exp(eta_Cl)
  pk1cptoral_(t=t, TIME=TIME, AMT=AMT, V=V, K=Cl/V, KA=ka)
})
# m1 <- RxODE::RxODE("
# ka = ka_pop * exp(eta_ka)
# V = V_pop * exp(eta_V)
# Cl = Cl_pop * exp(eta_Cl)
# d/dt(A0) = -ka*A0;
# d/dt(A1) = ka*A0 - Cl/V * A1;
# CONC = A1 / V;
# ")
theta <- MlxConnectors::getEstimatedPopulationParameters()
data <- MlxConnectors::getData()
data <- read_tsv(data$dataFile, na=".")

tdmorePred <- data %>% group_by(ID) %>% do({
  data <- .data %>% mutate(id=ID, time=TIME)
  regimen <- data %>% filter(!is.na(AMT)) %>% select(TIME, AMT)
  covariates <- theta[1:3]
  parameters <- data %>% left_join(mlxEta, by="id") %>%
    distinct(eta_ka, eta_V, eta_Cl) %>%
    unlist

  pred <- model_predict(m1,
                unlist(.$TIME),
                regimen=regimen,
                covariates=covariates,
                parameters = parameters)
  mutate(.data, IPRED=pred$CONC, time=TIME, id=ID)
})

test_that("Prediction using same EBE gives same results", {
  expect_equal(mlxPred %>% pull(IPRED),
               tdmorePred %>% filter(!is.na(CONC)) %>% pull(IPRED),
               tolerance=1E-12)
})

ggplot(mlxPred) +
  geom_line(aes(x=time, y=IPRED, color="Monolix")) +
  geom_line(data=tdmorePred, aes(x=time, y=IPRED, color="tdmore")) +
  facet_wrap(~id) +
  labs(title="Simulation comparison")
ggsave("mlx_theophyline_simul.png", width=16, height=9)


omega <- theta[4:6] ^ 2
names(omega) <- names(omega) %>% stringr::str_replace("omega", "eta")
m2 <- m1 %>% tdmore(
  #res_var=list( errorModel("CONC", add=theta['a'], prop=theta['b'], type="combined2") ),
  res_var=list( errorModel("CONC", add=theta['a'], prop=theta['b'], type="combined1") ),
  omega=omega,
  parameters=names(omega)
)

tdmorePred <- data %>% group_by(ID) %>% do({
  data <- .data %>% mutate(id=ID, time=TIME)
  regimen <- data %>% filter(!is.na(AMT)) %>% select(TIME, AMT)
  covariates <- theta[1:3]
  obs <- data %>% filter(!is.na(CONC)) %>% select(TIME, CONC)

  ipred <- estimate(m2, observed=obs, regimen=regimen, covariates=covariates)
  mlxParameters <- data %>% left_join(mlxEta, by="id") %>%
    distinct(eta_ka, eta_V, eta_Cl) %>%
    unlist

  tibble(
    ID=data$ID[1],
    ipred=list(ipred)
  )
}) %>% mutate(
  predicted=map(ipred, predict)
)

tdmoreEta <- tdmorePred %>% ungroup() %>% mutate(coef=map(ipred, ~as_tibble(t(coef(.x))))) %>% unnest(coef) %>%
  mutate(id=ID) %>% select(-ID, -ipred, -predicted) %>% select(!!colnames(mlxEta)) %>% as.data.frame

test_that("EBE estimation gives same estimates", {
  expect_equal(tdmoreEta, mlxEta, tolerance=1e-5)
})

tdmoreEta %>% gather(key, monolix, -id) %>% left_join(
  mlxEta %>% gather(key, tdmore, -id)
) %>% ggplot() +
  geom_point(aes(x=tdmore / monolix - 1, y=id)) +
  geom_vline(xintercept=0) +
  facet_wrap(~key) +
  labs(title="Relative difference between ETA estimates", subtitle="tdmore vs monolix")
ggsave("mlx_theophyline_eta.png", width=16, height=9)

test_that("EBE estimation gives same predictions", {
  expect_equal(mlxPred$IPRED,
               tdmorePred %>% unnest(predicted) %>% pull(CONC),
               tolerance=1e-6)
})

ggplot(mapping=aes(x=time, group=id)) +
  geom_point(data=mlxPred, aes(y=CONC)) +
  geom_line(data=mlxPred, aes(y=IPRED, color="Monolix")) +
  geom_line(data=tdmorePred %>% unnest(predicted) %>% mutate(id=ID), aes(x=TIME, y=CONC, group=ID, color="TDMore estim")) +
  facet_wrap(~id) +
  scale_y_log10()

ggplot(mapping=aes(x=id, y=value )) +
  geom_point(data=mlxEta %>% gather(key=eta, value=value, -id),
             aes(color="Monolix")) +
  geom_point(data= tdmoreEta %>% gather(key=eta, value=value, -id),
              aes(color="tdmore")) +
  geom_hline(yintercept=0) +
  facet_wrap(~eta)
