#' Create a TDM-capable model from an nlmixr UI object.
#'
#' @param model the nlmixr UI object
#' @param ... extra arguments will be passed to the model_predict call
#'
#' @return An object of class tdmore, which can be used to estimate posthoc Bayesian parameters.
#' The model contained within is actually an RxODE object.
#' @export
tdmore.nlmixrUI <- function(model, ...) {
  assert_that(class(model) %in% c("nlmixrUI"))

  # The processing below relies on the nlmixrUI object structure
  # It the structure changes, we'll have to adapt the code

  # Retrieving the parameters from omega matrix
  omega <- model$omega
  parameters <- colnames(omega)

  # Collecting all covariates
  covariates <- model$all.covs

  # Summarising the error models
  errorDf <- data.frame(cond=model$condition, errorType=model$err, value=model$est)
  errorDf <- errorDf %>% subset(errorType %in% c("add", "prop", "exp"))
  predDf <- model$predDf
  assert_that(nrow(predDf) >= 1, msg = "No error model defined, please define one")
  errorModels <- list()
  for (index in 1:nrow(predDf)) {
    row <- predDf[index,]
    add <- errorDf %>% subset(cond==as.character(row$cond) & errorType=="add")
    prop <- errorDf %>% subset(cond==as.character(row$cond) & errorType=="prop")
    exp <- errorDf %>% subset(cond==as.character(row$cond) & errorType=="exp")
    err <- errorModel(var = as.character(row$var),
                      add = if(nrow(add) > 0) as.numeric(add$value) else 0,
                      prop = if(nrow(prop) > 0) as.numeric(prop$value) else 0,
                      exp = if(nrow(exp) > 0) as.numeric(exp$value) else 0)
    errorModels[[length(errorModels) + 1]] <- err
  }

  # Model and structural equations processing
  modelFunction <- model$theta.pars
  thetaIndexes <- !is.na(model$ntheta)
  THETA <- model$est[thetaIndexes]

  modelCode <- as.character(body(modelFunction))
  modelCode <- modelCode[2:length(modelCode)] # Remove { character

  # Replace THETAS[x] variables by their real values
  index <- 1
  for(value in THETA) {
    modelCode <- gsub(paste0("THETA\\[", index, "\\]"), value, modelCode)
    index <- index + 1
  }

  # Remove all parameters assignments
  for(parameter in parameters) {
    modelCode <- modelCode[!startsWith(modelCode, paste0(parameter, " = "))]
  }

  # Merge model code with model$rxode code
  modelCode <- paste("   ", modelCode, collapse = "\n")
  rxOdeModelCode <- paste(modelCode, model$rxode, collapse = "\n")

  cat(rxOdeModelCode)

  # Create RxODE object
  rxModel <- RxODE::RxODE(rxOdeModelCode)

  # Construct the TDMore object
  tdmore <- structure(list(
    model=rxModel,
    omega=omega,
    res_var=errorModels,
    parameters=parameters,
    covariates=covariates,
    extraArguments=list(...)
  ), class="tdmore")

  # Check consistency and return
  return(checkTdmore(tdmore))
}
