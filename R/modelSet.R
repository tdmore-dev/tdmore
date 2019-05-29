#' Create a TDMore set.
#' Note that the order of the given models matters.
#' Order must be defined from the most restrictive model (with a high number of covariates) to the lowest restrictive model (with a low number of covariates).
#' These models will be tried one by one from left to right with the 'observed' data.
#' The first model matching the covariates from the 'observed' data will be used.
#'
#' @param ... 1 or more TDMore or TDMore mixture models, to be added in the set.
#'
#' @return a tdmore_set object
#' @export
tdmore_set <- function(...) {
  models <- list(...)
  for (model in models) {
    assert_that(is.tdmore(model) |
                  inherits(model, "tdmore_mixture"),
                msg = "Only tdmore or tdmore_mixture models can be added to a tdmore set")
  }
  assert_that(length(models) >= 1, msg = "You should provide at least one tdmore model")

  tdmoreSet <- structure(list(models = models), class = "tdmore_set")

  return(tdmoreSet)
}

#' Predict from a tdmore model set.
#'
#' @inheritParams predict.tdmore
#' @param ... extra arguments for the call to the model
#'
#'
#' @return A data.frame with all observed values at the given time points
#' @export
predict.tdmore_set <- function(object, newdata, regimen=NULL, parameters=NULL, covariates=NULL, se=FALSE, level=0.95, ...) {
  tdmore_set <- object
  chosenModel <- findFirstCompatibleModel(tdmore_set, covariates)
  return(chosenModel %>% predict(newdata, regimen, parameters, covariates, se, level, ...))
}

#' Find the first model in the TDMore set compatible with the specified covariates.
#'
#' @param tdmore_set a tdmore_set model
#' @param covariates the model covariates
#'
#' @return A tdmore or tdmore_mixture model
findFirstCompatibleModel <- function(tdmore_set, covariates) {
  chosenModel <- NULL
  for (model in tdmore_set$models) {
    isMixtureModel <- inherits(model, "tdmore_mixture")
    if (isMixtureModel) {
      # Take default tdmore model from the mixture
      modelToCheck <- model$models[[which.max(model$probs)[1]]]
    } else {
      # Model to check is a tdmore model
      modelToCheck <- model
    }
    if (!is.null(chosenModel)) break
    result = tryCatch({
      checkCovariates(modelToCheck, covariates)
      chosenModel <- model
    }, error = function(e) {})
  }
  assert_that(!is.null(chosenModel), msg = "No model is compatible with the provided covariates")
  return(chosenModel)
}
