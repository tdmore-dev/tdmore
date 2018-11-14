#' Create a TDM-capable model from a deSolve model
#'
#' @param parameters list of parameter names to be passed to the original function
#' @param add additive residual error, as stdev
#' @param prop proportional residual error, as stdev
#' @param exp exponential residual error, as stdev. The exponential error cannot be used in conjunction with the additive or proportional error
#' @param ... Arguments to deSolve
#'
#' @return An object of class tdmore, which can be used to estimate posthoc Bayesian parameters
#' @export
#'
#' @example inst/examples/deSolve.R
tdmore_deSolve <- function(parameters, add=0, prop=0, exp=0, ...) {
  assertthat::is.number(add)
  assertthat::is.number(prop)
  assertthat::is.number(exp)
  ## Exponential and add/prop are mutually exclusive
  if(exp != 0) assert_that(add==0 & prop==0)
  if(add != 0 || prop != 0) assert_that(exp == 0)

  structure(list(
    model=structure(list(extraArgs=list(...)), class="tdmore_deSolve"),
    res_var=list(add=add, prop=prop, exp=exp),
    parameters=parameters
  ), class="tdmore")
}

#' Predict for tdmore_deSolve
#'
#' @param model The model function
#' @param newdata data frame with at least a TIME column, and all observed data. The observed data will be compared to the model predictions.
#' If not specified, we estimate the population prediction
#' @param regimen dataframe with column 'TIME' and adhering to standard NONMEM specifications otherwise (columns AMT, RATE, CMT).
#' @param parameters either a dataframe with column 'TIME' and a column for each covariate and parameter, or a named numeric vector
#' @param covariates NOT IMPLEMENTED
#' @param extraArguments extra arguments to use
#'
#' @export
#'
#' @return
#' A data.frame similar to the observed data frame, but with predicted values.
#' @importFrom RxODE eventTable
#' @importFrom deSolve ode
#'
#' @keywords internal
model_predict.tdmore_deSolve <- function(model, newdata, regimen=data.frame(TIME=c()), parameters=c(), covariates=NULL, extraArguments=list()) {
  if(any(parameters > 10 | parameters < -10)) {
    warning("Solver requested with highly unlikely parameters, returning NA")
    newdata[, colnames(newdata) != "TIME"] <- NA
    return(newdata) #we cannot calculate this...
  }
  # Verify arguments are good
  samplingTimes <- NULL
  oNames <- NULL
  if("data.frame" %in% class(newdata)) {
    assert_that("data.frame" %in% class(newdata))
    assert_that("TIME" %in% colnames(newdata))
    oNames <- colnames(newdata)
    oNames <- oNames[oNames != "TIME"]
    assert_that(length(oNames) > 0)
    samplingTimes <- newdata$TIME
  } else {
    assert_that(is.numeric(newdata))
    samplingTimes <- newdata
  }

  assert_that("data.frame" %in% class(regimen))
  assert_that(all(c("TIME", "AMT") %in% colnames(regimen)))
  if("RATE" %in% colnames(regimen)) stop("Rate not allowed in deSolve model; use AMT and work this out yourself.")
  #assert_that(all(colnames(regimen) %in% c("TIME", "AMT", "RATE", "CMT", "II", "ADDL")))
  #if("II" %in% colnames(regimen) || "ADDL" %in% colnames(regimen))
  #  assert_that(all(c("II", "ADDL") %in% colnames(regimen)))
  assert_that(all(colnames(regimen) %in% c("TIME", "AMT", "CMT", "II", "ADDL")))

  params = NULL
  if(class(parameters) == "data.frame") {
    stop("Changing parameters is not supported in deSolve")
  } else {
    assert_that(is.numeric(parameters))
    pNames <- names(parameters)
    params = parameters
  }

  # All arguments look good, let's prepare the simulation
  regimen <- flatten(regimen)
  eventdat <- data.frame(
    var=1,
    time=regimen$TIME,
    value=regimen$AMT,
    method=2 #add
  )
  ## Events at time 0 are not applied in deSolve... *grmbl*
  #initevent <- subset(eventdat, time==0)
  #if(nrow(initevent) > 0) {
  #  y <- model$extraArgs$y
  #  for(i in seq_len(nrow(initevent))) {
  #    row <- initevent[i, ,drop=FALSE]
  #    y[row$var] <- y[row$var] + row$value
  #  }
  #  model$extraArgs$y <- y
  #}

  #eventdat <- subset(eventdat, time!=0)
  #eventdat <- subset(eventdat, min(samplingTimes) <= time & time <= max(samplingTimes))
  if(nrow(eventdat) == 0) eventdat = NULL

  result <- do.call(deSolve::ode, args=c(
    list(times=c(0, samplingTimes), parms=params, events=list(data=eventdat)),
    model$extraArgs
    ) )

  # Only get the values we want
  result <- as.data.frame(result)
  result <- result[ result$time %in% samplingTimes, ]
  result$TIME <- result$time
  if(!is.null(oNames)) result <- result[, c("TIME", oNames)]
  result
}
