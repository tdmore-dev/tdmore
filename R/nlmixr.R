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

  # Theta's and omega's names
  thetas <- names(model$theta)

  # Collecting all covariates
  covariates <- model$all.covs

  # Search for all parameters
  parameters <- model$rest.vars[!(model$rest.vars %in% thetas)]
  parameters <- parameters[!parameters %in% covariates]

  # Default values
  add <- 0
  prop <- 0
  exp <- 0

  # Collect the error model params
  addIndex <- which(model$err=="add")
  assert_that(length(addIndex) <= 1)
  if(length(addIndex) == 1) {
    add <- model$est[[addIndex]]
  }
  propIndex <- which(model$err=="prop")
  assert_that(length(propIndex) <= 1)
  if(length(propIndex) == 1) {
    prop <- model$est[[propIndex]]
  }
  expIndex <- which(model$err=="exp")
  assert_that(length(expIndex) <= 1)
  if(length(expIndex) == 1) {
    exp <- model$est[[expIndex]]
  }

  # Exponential and add/prop are mutually exclusive
  if(exp != 0) assert_that(add==0 & prop==0)
  if(add != 0 || prop != 0) assert_that(exp == 0)

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

  structure(list(
    model=rxModel,
    res_var=list(add=add, prop=prop, exp=exp),
    parameters=parameters,
    extraArguments=list(...)
  ), class="tdmore")
}
