## Check `r-source/src/library/stats/R/models.R` for a reference on all
## model functions in the R `stats` package. We try to use these function
## definitions as a base here.

assert_that <- assertthat::assert_that
# Structural model: how to predict --------------------------------------------------------
#' Prepare a cache object
#'
#' @param model The model itself.
#' @param times The times at which to generate predictions. May be empty, but never NULL.
#' @param regimen dataframe with column 'TIME' and adhering to standard NONMEM specifications otherwise (columns AMT, RATE, CMT)
#' @param parameters named vector
#' @param covariates named vector, or data.frame with column 'TIME', and at least TIME 0
#' @param iov IOV terms
#' @param extraArguments named list with extra arguments to use for call
#'
#' @return
#' A list of cached values
#'
#' @export
#' @keywords internal
#'
model_prepare <- function(model, times, regimen, parameters, covariates, iov, extraArguments) {
  UseMethod("model_prepare")
}

model_prepare.default <- function(model, times, regimen, parameters, covariates, iov, extraArguments) {
  NULL
}

#' Prototype predict function. Implement this to add your own model type to tdmore.
#'
#' @param model The model itself.
#' @param times The times at which to generate predictions. May be empty, but never NULL.
#' @param regimen dataframe with column 'TIME' and adhering to standard NONMEM specifications otherwise (columns AMT, RATE, CMT)
#' @param parameters named vector
#' @param covariates named vector, or data.frame with column 'TIME', and at least TIME 0
#' @param iov IOV terms
#' @param extraArguments named list with extra arguments to use for call
#'
#' @return
#' A data.frame with model predictions at each `times`
#'
#' @details
#' The model_predict function will generate a data.frame with the model output.
#' The data.frame will have a TIME column, and all model output in separate columns.
#'
#' In case the `times` vector is empty, it should generate an empty data.frame. The data.frame
#' should still contain all columns that would usually be generated.
#' Furthermore, it should check all arguments in the same way as usual, and even try to execute
#' a simulation.
#'
#' @export
#' @keywords internal
#'
model_predict <- function(model, times, regimen, parameters, covariates, iov, extraArguments, cache=NULL) {
  UseMethod("model_predict")
}

# Main API ----------------------------------------------------------------
#' tdmore is a generic function to append TDM functionality to a pharmacometrics model.
#'
#' In some cases, this function will automatically guess the required parameters from the model specified (e.g. names of parameters, names of covariates, residual error model).
#' In other cases, you will have to specify this manually.
#'
#' @param model the base model
#' @param ... extra arguments
#'
#' @return An object of class tdmore, which can be used to estimate posthoc Bayesian parameters
#' @export
tdmore <- function(model, ...) {
  UseMethod("tdmore")
}

#' Default function, catching unsupported types of models
#'
#' @param model model object
#' @param ... ignored
tdmore.default <- function(model, ...) {
  stop("Model of class ", class(model), " not supported.")
}

#' Test if object is of class `tdmore`
#'
#' @param x element to test
#' @export
is.tdmore <- function(x) {
  inherits(x, "tdmore")
}

#' Predict from a tdmore model
#'
#' @param object Object of class inheriting from `tdmore`
#' @param newdata Data.frame of new data with at least a TIME column, and additional columns.
#' Any columns that match the model output will be replaced by the predicted values.
#' Alternatively, a numeric vector used as TIME. The resulting data.frame will
#' contain all values predicted by the model.
#' @param regimen Treatment regimen, or NULL for no treatment
#' @param parameters A named numeric vector with the parameter values to use.
#' Any model parameters that are not provided are assumed `0`.
#' @param covariates a named numeric vector with the covariates to use.
#' Alternatively, a data.frame with column `TIME` and columns with time-varying covariates.
#' The first row of the data.frame should be time 0.
#' @inheritParams model.frame.tdmore
#' @param ... extra arguments for the call to the model
#'
#' @details
#' TODO: This is such an important function, that it merits a little more documentation
#'
#' @return A data.frame with all observed values at the given time points
#' @export
predict.tdmore <- function(object, newdata, regimen=NULL, parameters=NULL, covariates=NULL, se=FALSE, level=0.95, ...) {
  tdmore <- object
  checkCovariates(tdmore, covariates)

  # Process regimen
  if(is.null(regimen)) regimen <- data.frame(TIME=numeric(), AMT=numeric())

  # Process parameters
  par <- processParameters(parameters, tdmore, regimen)

  # Retrieve times vector from newdata
  if(is.data.frame(newdata)) times <- newdata$TIME
  else times <- as.numeric(newdata)

  # Call to model_predict
  predicted <- model_predict(model=tdmore$model, times=times, regimen=regimen, parameters=par, covariates=covariates, iov=tdmore$iov, extraArguments=c(..., tdmore$extraArguments), cache=tdmore$cache)

  if (is.data.frame(newdata)) {
    i <- colnames(newdata) %in% colnames(predicted)
    assert_that(all(i), msg=paste0("newdata contains unknown column(s): ",
                        paste(colnames(newdata)[!i], collapse="; ")))
    # Only use the outputs specified in newdata
    predicted <- predicted[ , colnames(newdata), drop=FALSE ]
  }

  # Finally, add residual error if needed
  result <- model.frame.tdmore(tdmore, predicted, se=se, level=level)
  result
}

#' Get the residual values of a predicted value vs the observed value
#'
#' @param object A TDMore object
#' @param observed The observed values, as data.frame
#' @param predicted What was predicted, as data.frame
#' @param weighted Should the residuals be divided by the standard deviation of the error model?
#' @param ... ignored
#'
#' @return A numeric vector, equivalent to IRES or IWRES
#' @export
#'
#' @importFrom stats dnorm
residuals.tdmore <- function(object, observed, predicted, weighted=FALSE, ...) {
  tdmore <- object
  result <- c()
  for (err in tdmore$res_var) {
    var <- err$var

    ipred <- predicted[, var, drop=TRUE]
    obs <- observed[, var, drop=TRUE]

    if(weighted) {
      tmp <- err$wres(ipred, obs)
    } else {
      tmp <- err$res(ipred, obs)
    }
    result <- c(result, tmp)
  }
  return(result)
}

#' Get the specified data values, with the upper and lower bounds provided by the residual error model
#'
#' @param formula a TDMore object
#' @param data data.frame with at least a TIME column
#' Or alternatively, a numeric vector.
#' Or NULL, to get an empty data.frame with TIME and a column for every residual error variable.
#' @param se TRUE to generate an additional xx.upper and xx.lower column for every xx column in the data dataset
#' @param level numeric value to specify the confidence interval
#' @param ... ignored
#'
#' @importFrom stats qnorm
#'
#' @return
#' a data.frame similar to data. It contains at least column TIME.
#' If a numeric vector was specified as `data`, the output is a data.frame with TIME column and additional NA-filled columns corresponding to the residual error variables.
#' If `se` is TRUE, the data.frame has additional columns xx.lower and xx.upper for all columns that match the residual error model.
#' @export
model.frame.tdmore <- function(formula, data, se=FALSE, level=0.95, ...) {
  tdmore <- formula
  if(is.null(data)) data <- numeric()

  if(is.data.frame(data)) {
    assert_that("TIME" %in% colnames(data))
  } else {
    # Built a data.frame with TIME column and all residual error columns
    data <- as.numeric(data)
    data <- data.frame(TIME=data)
    for(err in tdmore$res_var) data[, err$var] <- numeric()
  }

  if(!se) return(data)

  oNames <- colnames(data)
  oNames <- oNames[oNames != "TIME"]

  if(is.na(level)) {
    ## Simply add RE
    for (err in tdmore$res_var) {
      var <- err$var
      if (!(var %in% oNames)) next
      obs <- data[, var, drop=TRUE]
      sd <- err$sigma(obs)
      q <- stats::rnorm(1)
      data[, var] <- obs + sd*q
    }
    return(data)
  }

  ## Calculate limits
  a <- (1 - level) / 2
  q <- qnorm(a) * (-1) # To have a positive value

  for (err in tdmore$res_var) {
    var <- err$var
    if (!(var %in% oNames)) next
    obs <- data[, var, drop=TRUE]
    sd <- err$sigma(obs)
    data[, paste0(var, ".lower")] <- obs - sd*q
    data[, paste0(var, ".upper")] <- obs + sd*q
  }
  return(data)
}


#' Print a tdmore object.
#'
#' @param x a tdmore object
#' @param ... ignored
#'
#' @export
print.tdmore <- function(x, ...) {
  cat("Structural model:", class(x$model), "\n")
  cat("Parameters:", x$parameters, "\n")
  cat("Covariates:", if(length(x$covariates)==0) "/" else x$covariates, "\n")
  cat("Output(s):\n")
  for (err in x$res_var) {
    print(err)
  }
  invisible(x)
}

#' Summarise a tdmore object.
#'
#' @param object a tdmore object
#' @param ... ignored
#'
#' @export
summary.tdmore <- function(object, ...) {
  structure(list(
    tdmore=object
  ), class="summary.tdmore")
}

#' Print a tdmore summary
#'
#' @param x a tdmore summary
#' @param ... ignored
#'
#' @export
print.summary.tdmore <- function(x, ...) {
  x <- x$tdmore
  cat("Structural model:", class(x$model), "\n\n")
  parameters <- data.frame(name=x$parameters, var=diag(x$omega), cv=sqrt(diag(x$omega)))
  cat("Parameters:\n")
  print(parameters, row.names = FALSE)
  cat("\n")
  cat("Covariates:", if(length(x$covariates)==0) "/" else x$covariates, "\n\n")
  errorDf <- data.frame(name=character(0), additiveError=numeric(0), proportionalError=numeric(0), exponentialError=numeric(0))
  for (index in 1:length(x$res_var)) {
    err <- x$res_var[[index]]
    errDf <- data.frame(name=err$var, type=err$type, additiveError=err$add, proportionalError=err$prop)
    errorDf <- rbind(errorDf, errDf)
  }
  cat("Residual error model:\n")
  print(errorDf, row.names = FALSE)
}

#' Check input parameters and return a standardised named num array with the parameter names and values.
#' If the initial values are not provided, zeroes will be used (=population values for ETAS)
#' Provided parameters will overwrite initial values.
#'
#' @param parameters user inputted array of parameters
#' @param tdmore the tdmore object
#' @param regimen the regimen
#' @param defaultValues the default values, in the right order
#'
#' @return a standardised named num array with the parameter names and values
#'
processParameters <- function(parameters, tdmore, regimen, defaultValues=NULL) {
  parameterNames <- getParameterNames(tdmore, regimen)

  if(is.null(defaultValues)) {
    defaultValues <- rep(0, length(parameterNames)) # start with population
    names(defaultValues) <- parameterNames
  }

  ## Use the defaultValues, but possibly extend with 0 for e.g. IOV
  assert_that(identical(
    names(defaultValues),
    parameterNames[seq_along(defaultValues)] ))
  Nmissing <- length(parameterNames) - length(defaultValues)
  par <- c(defaultValues, rep(0, Nmissing))
  names(par) <- parameterNames

  if(is.null(parameters)) return(par)
  if(length(parameters)==1 && is.na(parameters)) return(par)

  assert_that(is.numeric(parameters))
  assert_that(all(names(parameters) %in% names(par)),
              msg=paste("Unknown parameters", paste(names(parameters[!(names(parameters) %in% names(par))]), collapse = ",")))
  par[names(parameters)] <- parameters # set par from argument
  iov <- tdmore$iov
  if(!is.null(iov)) {
    for(iov_term in iov) {
      updatedPar <- parameters[which(names(parameters)==iov_term)]
      if(length(updatedPar) > 0) {
        parIndexes <- which(names(par)==iov_term)
        # assert_that(length(updatedPar) == length(parIndexes),
        #             msg=paste0("Incorrect number of initial values for IOV term ", iov_term, " (", length(parIndexes), " needed)"))
        par[parIndexes[1:length(updatedPar)]] <- updatedPar
      }
    }
  }
  return(par)
}

#' Get the parameter names. If IOV is defined in the model,
#' IOV terms will be added at the end of the character array and repeated as many times as there are occasions in the regimen.
#' E.g. 2 occasions: "ECL", "EKA", "ECL_IOV", "EKA_IOV", "ECL_IOV", "EKA_IOV"
#'
#' @param tdmore the tdmore object
#' @param regimen the regimen
#' @return an ordered character array with the parameters names
getParameterNames <- function(tdmore, regimen) {
  iov <- tdmore$iov # IOV terms order is same as in parameters
  parameters <- tdmore$parameters
  if(is.null(iov)) {
    retValue <- tdmore$parameters
  } else {
    retValue <- c(parameters[!(parameters %in% iov)],
                  rep(parameters[(parameters %in% iov)], getMaxOccasion(regimen)))
  }
  return(retValue)
}
