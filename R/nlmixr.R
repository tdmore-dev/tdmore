#' @param iov list of parameter names related to IOV, NULL if no IOV
#' @name tdmore
#' @export
tdmore.nlmixrUI <- function(model, iov=NULL, ...) {
  stopifnot(inherits(model, "nlmixrUI"))

  # The processing below relies on the nlmixrUI object structure
  # It the structure changes, we'll have to adapt the code

  # Retrieving the parameters from omega matrix
  omega <- model$omega
  parameters <- colnames(omega)

  # Collecting all covariates
  covariates <- model$all.covs

  # Summarising the error models
  errorDf <- data.frame(cond=model$condition, errorType=model$err, value=model$est)
  errorDf <- errorDf %>% subset(errorDf$errorType %in% c("add", "prop", "exp"))
  predDf <- model$predDf
  if(nrow(predDf) == 0) stop("No error model defined, please define one")
  errorModels <- list()
  for (index in 1:nrow(predDf)) {
    row <- predDf[index,]
    add <- errorDf %>% subset(errorDf$cond==as.character(row$cond) & errorDf$errorType=="add")
    prop <- errorDf %>% subset(errorDf$cond==as.character(row$cond) & errorDf$errorType=="prop")
    err <- errorModel(var = as.character(row$var),
                      add = if(nrow(add) > 0) as.numeric(add$value) else 0,
                      prop = if(nrow(prop) > 0) as.numeric(prop$value) else 0)
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
  rxodeCode <- gsub(paste0("\\n(cmt|dvid)\\(.*\\);\\n"), "", model$rxode)
  rxOdeModelCode <- paste(modelCode, rxodeCode, collapse = "\n")

  # Create RxODE object
  rxModel <- RxODE::RxODE(rxOdeModelCode)

  # Construct the TDMore object
  tdmore <- structure(list(
    model=rxModel,
    omega=omega,
    res_var=errorModels,
    parameters=parameters,
    covariates=covariates,
    iov=iov,
    extraArguments=list(...)
  ), class="tdmore")

  # Check consistency and return
  return(checkTdmore(tdmore))
}

#' @name tdmore
#' @export
tdmore.nlmixrFitCore <- function(model, ...) {
  # Use the included nlmixrUI object
  tdmore(model$uif)
}
