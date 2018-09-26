assert_that <- assertthat::assert_that
# Structural model: how to predict --------------------------------------------------------
#' Prototype predict function
#'
#' @param model The model itself.
#' @param newdata dataframe with at least a column 'TIME' and other values. A prediction will be generated for each value.
#' @param regimen dataframe with column 'TIME' and adhering to standard NONMEM specifications otherwise (columns AMT, RATE, CMT)
#' @param parameters named vector
#' @param covariates named vector, or data.frame with column 'TIME', and at least TIME 0
#' @param extraArguments named list with extra arguments to use for call
#'
#' @return
#' A data.frame similar to the observed data frame, but with predicted values.
#' @export
#'
model_predict <- function(model, newdata, regimen, parameters, covariates, extraArguments) {
  UseMethod("model_predict")
}

# Main API ----------------------------------------------------------------
#' Create a new TDM model
#'
#' @param model The model itself
#' @param parameters character vector of parameter names, or NULL to auto-detect
#' @param add additive residual error, as stdev
#' @param prop proportional residual error, as stdev
#' @param exp exponential residual error, as stdev. The exponential error cannot be used in conjunction with the additive or proportional error
#' @param ... extra arguments
#'
#' @return An object of class tdmore, which can be used to estimate posthoc Bayesian parameters
#' @export
#'
#' @example inst/examples/RxODE.R
tdmore <- function(model, parameters=NULL, add=0, prop=0, exp=0, ...) {
  UseMethod("tdmore")
}

#' Default function, catching unsupported types of models
#'
#' @param model model object
#' @param ... ignored
tdmore.default <- function(model, ...) {
  stop("Model of class ", class(model), " not supported.")
}

#' Predict from a tdmore model
#'
#' @param object TDMore object
#' @param newdata Data.frame of new data with at least a TIME column, and blank columns with all desired predictions
#' or a numeric vector to predict all possible values from the model
#' @param regimen Treatment regimen
#' @param parameters The parameter values to use
#' missing values are taken from the population
#' @param se Whether to add residual error
#' @param level How much residual error to add
#' @param ... ignored
#'
#' @seealso tdmore::model_predict
#'
#' @return A data.frame with all observed values at the given time points
#' @export
predict.tdmore <- function(object, newdata, regimen, parameters=NULL, covariates=NULL, se=FALSE, level=0.95, ...) {
  tdmore <- object
  pars <- rep(0, length(tdmore$parameters)) #population prediction
  names(pars) <- tdmore$parameters
  if(!is.null(parameters)) {
    assert_that(is.numeric(parameters))
    assert_that(all(names(parameters) %in% names(pars)))
    pars[names(parameters)] <- parameters ## set pars from argument
  }

  predicted <- model_predict(model=tdmore$model, newdata=newdata, regimen=regimen, parameters=pars, covariates=covariates, extraArguments=tdmore$extraArguments)
  assert_that("data.frame" %in% class(predicted))
  result <- model.frame.tdmore(tdmore, predicted, se=se, level=level)
  result
}

#' Get the residual values of a predicted value vs the observed value
#' This is based on the log-likelihood of the observed value in the residual error model
#'
#' @param object A TDMore object
#' @param observed The observed values, as data.frame
#' @param predicted What was predicted, as data.frame
#' @param log if TRUE, report as logPdf
#' @param ... ignored
#'
#' @return A numeric vector
#' @export
#'
#' @importFrom stats dnorm
residuals.tdmore <- function(object, observed, predicted, log=TRUE, ...) {
  tdmore <- object
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

#' Get the specified observed values, with the upper and lower bounds provided by the residual error model
#'
#' @param formula a TDMore object
#' @param observed data.frame with at least a TIME column
#' @param se TRUE to generate an additional xx.upper and xx.lower column for every xx column in the observed dataset
#' @param level numeric value to specify the confidence interval
#' @param ... ignored
#'
#' @importFrom stats qnorm
#'
#' @return a data.frame similar to observed, but with additional columns xx.lower and xx.upper (in case of se=TRUE)
#' @export
model.frame.tdmore <- function(formula, observed, se=FALSE, level=0.95, ...) {
  tdmore <- formula
  if(!se) return(observed)
  if(is.null(observed)) return(NULL)
  assert_that("data.frame" %in% class(observed))
  assert_that("TIME" %in% colnames(observed))
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

#' Get the original model
#'
#' @param x TDMore object
#' @param ... ignored
#'
#' @return The original model used when defining the TDMore object
#' @export
formula.tdmore <- function(x, ...) {x$model}
