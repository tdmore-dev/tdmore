#' Create a TDMore set.
#'
#' @param ... 1 or more tdmore models, to be added in the set.
#'
#' @return a tdmore_set object
#' @export
tdmore_set <- function(...) {
  models <- list(...)
  for (model in models) {
    assert_that("tdmore" %in% class(model), msg = "Only tdmore models can be added to the set")
  }
  assert_that(length(models) >= 1, msg = "You should provide at least one tdmore model")

  tdmoreSet <- structure(list(
    models=models
  ), class="tdmore_set")

  # Check consistency and return
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
#' @return A data.frame with all observed values at the given time points
#' @export
findFirstCompatibleModel <- function(tdmore_set, covariates) {
  chosenModel <- NULL
  for (model in tdmore_set$models) {
    if (!is.null(chosenModel)) break
    result = tryCatch({
      checkCovariates(model, covariates)
      chosenModel <- model
    }, error = function(e) {})
  }
  assert_that(!is.null(chosenModel), msg = "No model is compatible with the provided covariates")
  return(chosenModel)
}


