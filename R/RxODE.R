#' Create a TDM-capable model from an RxODE model
#'
#' @param model the RxODE model
#' @param res_var the residual variability
#' @param parameters list of parameter names, or NULL to use all parameters names from RxODE
#' @param omega omega variance-covariance matrix, or NULL to use a diagonal matrix of variances 1
#' @param iov list of parameter names related to IOV, NULL if no IOV
#' @param ... extra arguments will be passed to the model_predict call
#'
#' @return An object of class tdmore, which can be used to estimate posthoc Bayesian parameters
#' @export
#'
#' @example inst/examples/RxODE.R
tdmore.RxODE <- function(model, res_var, parameters=NULL, omega=NULL, iov=NULL, ...) {
  tdmore <- structure(list(
    model=model,
    omega=omega,
    res_var=res_var,
    parameters=parameters,
    covariates=NULL, # Computed automatically in checkTdmore
    iov=iov,
    extraArguments=list(...)
  ), class="tdmore")

  # Check consistency and return
  # TODO: checkTdmore does more than just check! It modifies!
  return(checkTdmore(tdmore))
}

#' Predict for RxODE
#'
#' @param model The model itself.
#' @param newdata
#' A dataframe with at least a column 'TIME' and other values that should be observed,
#' or a numeric TIME vector to produce all values that can be observed by this model
#' @param regimen dataframe with column 'TIME' and adhering to standard NONMEM specifications otherwise (columns AMT, RATE, CMT).
#' @param parameters named numeric vector with the parameters
#' @param covariates named vector, or data.frame with column 'TIME', and at least TIME 0
#' @param iov character array with the IOV terms, can be NULL
#' @param extraArguments extra arguments to use
#'
#' @export
#'
#' @return
#' A data.frame similar to the observed data frame, but with predicted values.
#' @importFrom RxODE eventTable rxSolve
#' @importFrom dplyr transmute left_join
#'
#' @keywords internal
#'
model_predict.RxODE <- function(model, times, regimen=data.frame(TIME=numeric()), parameters=numeric(), covariates=NULL, iov=NULL, extraArguments=list()) {
  ### RxODE sometimes errors out...
  ### Probably not a solver issue, but rather
  ### issue that the DLL of RxODE is no longer loaded
  ### (especially when running in Simulo)
  if(!RxODE::rxDllLoaded(model)) {
    warning("RxDLL was not loaded, loading...")
    RxODE::rxLoad(model)
    if(!RxODE::rxDllLoaded(model)) stop("Tried to reload RxDLL, but failed. Cannot continue!")
  }
  modVars <- model$get.modelVars()

  # Check times and regimen objects
  checkTimes(times)
  checkRegimen(regimen, iov)

  # IOV processing
  iovPrediction <- !is.null(iov)
  iovIndexes <- which(names(parameters) %in% iov)
  iovParameters <- parameters[iovIndexes]
  parameters <- parameters[!duplicated(names(parameters))]

  # Flatten the regimen (additional doses are converted)
  # This is done only for IOV, because RxODE is slower with a flattened regimen
  if(iovPrediction) regimen <- flatten(regimen)

  # Set up the parameters
  params = NULL
  assert_that(is.numeric(parameters))
  pNames <- names(parameters) # parameter names
  params = parameters # parameters + fixed covariate values
  cNames <- c()

  covariateAsDataFrame <- "data.frame" %in% class(covariates)
  if(covariateAsDataFrame) {
    assert_that("TIME" %in% colnames(covariates))
    assert_that(0 %in% covariates$TIME)
    cNames <- colnames(covariates)
    cNames <- cNames[cNames != "TIME"]
  } else if (is.numeric(covariates)){
    cNames <- names(covariates)
    params = c(parameters, covariates)
    covs = NULL
  } else if (length(covariates) == 0) {
    covs=NULL
  } else {
    stop("Covariates in wrong format")
  }

  pNames <- c(pNames, cNames) # parameter names + all covariate names
  assert_that(all(pNames %in% modVars$params))
  assert_that(all(modVars$params %in% pNames))

  # Occasion processing
  occasionTimes <- getOccasionTimes(regimen)
  occasions <- getMaxOccasion(regimen)

  # ODE's initial values
  inits <- rep(0, length(modVars$state))
  names(inits) <- modVars$state
  retValue <- NULL

  # RxODE does not allow to simulate 'nothing'
  # Manually construct an empty data.frame with the right columns
  if(length(times) == 0 && nrow(regimen) == 0) {
    result <- data.frame()
    for(i in c("TIME", modVars$lhs, modVars$state)) result[, i] <- numeric()
    return(result)
  }

  # Predict each occasion
  for(occasion in 1:occasions) {
    last <- occasion == occasions
    currentCovariates <- covariates

    if(iovPrediction) {
      occasionTime <- occasionTimes[occasion]
      nextOccasionTime <- if(last) {Inf} else {occasionTimes[occasion + 1]}
      currentRegimen <- regimen %>% subset(regimen$TIME >= occasionTime & regimen$TIME < nextOccasionTime)
      currentTimes <- times[times >= occasionTime & times <= nextOccasionTime]
      currentTimes <- if(is.finite(nextOccasionTime)) {unique(c(currentTimes, nextOccasionTime))} else {currentTimes}
      if(covariateAsDataFrame) {
        currentCovariates <- currentCovariates %>% subset(currentCovariates$TIME >= occasionTime & currentCovariates$TIME < nextOccasionTime)
      }
      for(iov_term in iov) {
        iovValueIndexes <- which(names(iovParameters)==iov_term)
        iovValue <- iovParameters[[iovValueIndexes[occasion]]]
        if(is.null(iovValue)) {
          stop(paste("Missing IOV values for", iov_term))
        }
        params[iov_term] <- iovValue
      }
    } else {
      currentRegimen <- regimen
      currentTimes <- times
    }

    # All arguments look good, let's prepare the simulation
    ev <- RxODE::eventTable()
    if(length(currentTimes) > 0) ev$add.sampling(time=currentTimes)
    ev <- addRegimenToEventTable(ev, currentRegimen)

    # Treat missing times
    if(covariateAsDataFrame) {
      missingTimes <- currentCovariates$TIME[ !(currentCovariates$TIME %in% ev$get.sampling()$time) ]
      if(length(missingTimes) > 0) ev$add.sampling(missingTimes)

      eventTable <- ev$get.EventTable()
      covs <-  eventTable %>% dplyr::transmute(TIME=eventTable$time) %>% dplyr::left_join(currentCovariates, by="TIME")
      for(i in colnames(covs))
        covs[, i] <- zoo::na.locf(covs[,i]) # last observation carried forward for NA values
    }

    # Run the simulation
    result <- do.call(RxODE::rxSolve, c( list(object=model, events=ev, params=params, covs=covs, inits=inits), extraArguments))

    # Adapt initial values for next iteration
    inits <- as.numeric(result[nrow(result),][modVars$state])
    names(inits) <- modVars$state

    # Remove last row if not last iteration
    if(!last) result <- result[-nrow(result),]

    # Only get the values we want
    result <- as.data.frame(result, row.names=NULL)
    names(result)[names(result)=="time"] <- "TIME"

    # Remove spurious sampling times due to covariates
    # But keep this a data.frame!
    result <- subset( result, result$TIME %in% times, drop=FALSE)

    # Bind to global result
    retValue <- rbind(retValue, result)
  }

  retValue
}

#' Add the given regimen into the RxODE event table.
#'
#' @param eventTable the RxODE event table
#' @param regimen the specified regimen
#'
#' @return a completed event table
#'
addRegimenToEventTable <- function(eventTable, regimen) {
  for(i in seq_len(nrow(regimen))) {
    row <- regimen[i, ,drop=FALSE]
    dosing.to = 1
    rate = NULL
    nbr.doses = 1
    dosing.interval = 24

    # dosing interval not used since regimen has been flattened
    if(isTRUE(row$II > 0)) {
      dosing.interval <- row$II
      nbr.doses <- row$ADDL
    }

    if(isTRUE(is.finite(row$RATE))) rate=row$RATE
    if(isTRUE(is.finite(row$DURATION))) {
      if(!is.null(rate)) stop("Cannot specify RATE and DURATION in the same treatment row ", i)
      rate = row$AMT / row$DURATION
    }

    if("CMT" %in% names(row)) dosing.to <- row$CMT
    eventTable$add.dosing(start.time = row$TIME,
                  dose=row$AMT,
                  dosing.to=dosing.to,
                  rate=rate,
                  dosing.interval=dosing.interval,
                  nbr.doses=nbr.doses)
  }
  return(eventTable)
}
