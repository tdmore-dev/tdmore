#' Create a TDMore mixture model.
#'
#' @param ... 2 or more tdmore models that describe different subpopulations
#' @param probs 'a priori' probabilities for belonging to the different subpopulations, numeric vector
#'
#' @return a tdmore_set object
#' @export
tdmore_mixture <- function(..., probs) {
  models <- list(...)
  for (model in models) {
    assert_that(inherits(model, "tdmore"), msg = "Only tdmore models can be added to the tdmore mixture")
  }
  assert_that(length(models) >= 2, msg = "You should provide at least two tdmore models to create a tdmore mixture")
  assert_that(is.numeric(probs), msg = "probs is not numeric")
  assert_that(length(models) == length(probs), msg = "There must be as many models as probabilities defined in the numeric vector 'probs'")
  assert_that(isTRUE(all.equal(target = 1, current = sum(probs), tolerance = 1e-6)), msg = "Sum of probabilities must be 1")
  tdmoreMixture <- structure(list(
    models=models,
    probs=probs
  ), class="tdmore_mixture")

  return(tdmoreMixture)
}

#' Predict from a TDMore mixture model. The model with the highest probability is used for the predictions.
#'
#' @inheritParams predict.tdmore
#' @param ... extra arguments for the call to the model
#'
#'
#' @return A data.frame with all observed values at the given time points
#' @export
predict.tdmore_mixture <- function(object, newdata, regimen=NULL, parameters=NULL, covariates=NULL, se=FALSE, level=0.95, ...) {
  tdmoreMixture <- object
  chosenModel <- tdmoreMixture$models[[which.max(tdmoreMixture$probs)]]
  return(chosenModel %>% predict(newdata, regimen, parameters, covariates, se, level, ...))
}

#' Calculate the Empirical Bayesian Estimate parameters that predict
#' a given dataset using a mixture model.
#'
#' @inheritParams estimate
#' @param ... extra parameters to pass to the optim function
#'
#' @return A tdmorefit object
#' @importFrom stats optim
estimateMixtureModel <- function(object, observed=NULL, regimen, covariates=NULL, par=NULL, method="L-BFGS-B", lower=NULL, upper=NULL, ...) {
  mixture <- object
  models <- mixture$models
  probs <- mixture$probs
  fits <- lapply(models, FUN=function(model) {
    estimate(model, observed=observed, regimen=regimen, covariates=covariates, par=par, method=method, lower=lower, upper=upper, ...)
  })
  mixtureProbs <-
    data.frame(lik = as.numeric(lapply(fits, FUN = function(fit) {exp(fit$logLik)})),
               probs = probs)
  mixtureProbs$IPkNumerator <- mixtureProbs$lik * mixtureProbs$probs
  mixtureProbs$IPk <- mixtureProbs$IPkNumerator/sum(mixtureProbs$IPkNumerator)
  winnerIndex <- which(mixtureProbs$IPk == max(mixtureProbs$IPk))

  return(fits[[winnerIndex]])
}

