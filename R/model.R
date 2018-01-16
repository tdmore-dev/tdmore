assert_that <- assertthat::assert_that
# Structural model: how to predict --------------------------------------------------------
#' Prototype predict function
#'
#' @param model The model itself.
#' @param observed dataframe with at least a column 'TIME' and other values. A prediction will be generated for each filled-in value.
#' @param regimen dataframe with column 'TIME' and adhering to standard NONMEM specifications otherwise (columns AMT, RATE, CMT)
#' @param parameters dataframe with column 'TIME' and a column for each covariate and parameter
#'
#' @return
#' A data.frame similar to the observed data frame, but with predicted values.
#'
model_predict <- function(model, observed, regimen, parameters) {
  UseMethod("model_predict")
}

#' Predict for RxODE
#'
#' @param model The model itself.
#' @param observed dataframe with at least a column 'TIME' and other values. A prediction will be generated for each filled-in value.
#' @param regimen dataframe with column 'TIME' and adhering to standard NONMEM specifications otherwise (columns AMT, RATE, CMT)
#' @param parameters either a dataframe with column 'TIME' and a column for each covariate and parameter, or a named numeric vector
#'
#' @return
#' A data.frame similar to the observed data frame, but with predicted values.
#'
model_predict.RxODE <- function(model, observed, regimen=data.frame(TIME=c()), parameters=c()) {
  modVars <- model$get.modelVars()

  # Verify arguments are good
  assert_that("data.frame" %in% class(observed))
  assert_that("TIME" %in% colnames(observed))
  oNames <- colnames(observed)
  oNames <- oNames[oNames != "TIME"]
  assert_that(length(oNames) > 0)
  assert_that(all(oNames %in% c(modVars$lhs, modVars$state)))

  assert_that("data.frame" %in% class(regimen))
  assert_that(all(c("TIME", "AMT") %in% colnames(regimen)))
  assert_that(all(colnames(regimen) %in% c("TIME", "AMT", "RATE", "CMT")))

  covs = NULL
  params = NULL
  if(class(parameters) == "data.frame") {
    assert_that("TIME" %in% colnames(parameters))
    pNames <- colnames(parameters)
    pNames <- pNames[pNames != "TIME"]
    covs = parameters
  } else {
    assert_that(is.numeric(parameters))
    pNames <- names(parameters)
    params = parameters
  }
  assert_that(all(pNames %in% modVars$params))
  assert_that(all(modVars$params %in% pNames))

  # All arguments look good, let's prepare the simulation
  ev <- RxODE::eventTable()
  ev$add.sampling(time=observed$TIME)
  for(i in 1:nrow(regimen)) {
    row <- regimen[i, ,drop=FALSE]
    dosing.to = 1
    rate = NULL
    if("RATE" %in% names(row)) rate = row$RATE
    if("CMT" %in% names(row)) dosing.to <- row$CMT
    ev$add.dosing(start.time = row$TIME,
                  dose=row$AMT,
                  dosing.to=dosing.to,
                  rate=rate)
  }

  # Run the simulation
  result <- model$solve(events=ev, params=params, covs=covs)

  # Only get the values we want
  result <- as.data.frame(result)
  result$TIME <- result$time
  result[, colnames(observed)]
}


# Main API ----------------------------------------------------------------
#' Create a new TDM model
#'
#' @param model
#' @param parameters character vector of parameter names, or NULL to auto-detect
#'
#' @return
#' a tdmore object
#' @export
#'
#' @examples
tdmore <- function(model, parameters=NULL, add=0, prop=0, exp=0) {
  assert_that(class(model) %in% c("RxODE")) #currently, only RxODE is supported

  # Check that parameters + covariates together supplies the parameters needed for the model
  modVars <- model$get.modelVars()
  if(is.null(parameters)) parameters=modVars$params
  assert_that( all.equal( sort(parameters), sort(modVars$params)) )

#  if(is.null(covariates)) covariates = modVars$params[ ! modVars$params %in% parameters ] #the rest
#  assert_that( all.equal( sort(c(parameters, covariates)), sort(modVars$params)) )

  assertthat::is.number(add)
  assertthat::is.number(prop)
  assertthat::is.number(exp)
  ## Exponential and add/prop are mutually exclusive
  if(exp != 0) assert_that(add==0 & prop==0)
  if(add != 0 || prop != 0) assert_that(exp == 0)

  structure(list(
    model=model,
    res_var=list(add=add, prop=prop, exp=exp),
    parameters=parameters
  ), class="tdmore")
}


predict.tdmore <- function(tdmore, observed, regimen, parameters, se=FALSE, level=0.95) {
  predicted <- model_predict(tdmore$model, observed, regimen, parameters)
  model.frame.tdmore(tdmore, predicted, se=se, level=level)
}

residuals.tdmore <- function(tdmore, observed, predicted, log=TRUE) {
  oNames <- colnames(observed)
  oNames <- oNames[oNames != "TIME"]
  i <- oNames %in% colnames(predicted)
  i <- oNames[i]
  ipred <- predicted[,i]
  obs <- observed[,i]
  res <- tdmore$res_var
  if(res$exp != 0) {
    sd <- res$exp
    dnorm(log(ipred), log(obs), sd, log=log) #calculate residual in the log domain
  } else {
    sd <- sqrt(res$add**2 + (res$prop * ipred)**2)
    ifelse(sd == 0 & ipred==obs, #treat points with 0 res. error special
           0, #ignore point
           dnorm(ipred, obs, sd, log=log) # calculate residual
    )
  }
}

model.frame.tdmore <- function(tdmore, observed, se=FALSE, level=0.95) {
  if(!se) return(observed)
  oNames <- colnames(observed)
  oNames <- oNames[oNames != "TIME"]
  a <- (1 - level) / 2
  q <- qnorm(a)
  res <- tdmore$res_var
  for(n in oNames) {
    obs <- observed[, n]
    if(res$exp != 0) {
      observed[, paste0(n, ".upper")] <- obs * exp(res$exp * q)
      observed[, paste0(n, ".lower")] <- obs * exp(-res$exp * q)
    } else {
      observed[, paste0(n, ".upper")] <- obs * (1 + res$prop*q) + res$add*q
      observed[, paste0(n, ".lower")] <- obs * (1 - res$prop*q) - res$add*q
    }
  }

  observed
}

formula.tdmore <- function(tdmore) {tdmore$model}
