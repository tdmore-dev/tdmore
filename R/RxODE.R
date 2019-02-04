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
#' @param iov character array with the parameters that should change for every different OCC, can be NULL
#' @param extraArguments extra arguments to use
#'
#' @export
#'
#' @return
#' A data.frame similar to the observed data frame, but with predicted values.
#' @importFrom RxODE eventTable rxSolve
#' @importFrom dplyr transmute left_join bind_rows arrange mutate_all
#' @importFrom rlang .data
#'
#' @keywords internal
#'
model_predict.RxODE <- function(model, times, regimen=data.frame(TIME=numeric()), parameters=numeric(), covariates=NULL, iov=NULL, extraArguments=list()) {
  modVars <- model$get.modelVars()
  # RxODE does not allow to simulate 'nothing'
  # Manually construct an empty data.frame with the right columns
  if(length(times) == 0) {
    result <- data.frame()
    for(i in c("TIME", modVars$lhs, modVars$state)) result[, i] <- numeric()
    return(result)
  }

  ### RxODE sometimes errors out...
  ### Probably not a solver issue, but rather
  ### issue that the DLL of RxODE is no longer loaded
  ### (especially when running in Simulo)
  if(!RxODE::rxDllLoaded(model)) {
    warning("RxDLL was not loaded, loading...")
    RxODE::rxLoad(model)
    if(!RxODE::rxDllLoaded(model)) stop("Tried to reload RxDLL, but failed. Cannot continue!")
  }

  # Check times and regimen objects
  checkTimes(times)
  checkRegimen(regimen, iov)

  # Check parameters
  assert_that(is.numeric(parameters))
  i <- names(parameters) %in% iov
  iovParameters <- parameters[i] #split the parameters in two arrays
  parameters <- parameters[!i]

  # Set up the covariates
  if (is.numeric(covariates)){
    parameters <- c(parameters, covariates)
    covariates <- NULL
  } else if(is.data.frame(covariates)) {
    assert_that("TIME" %in% colnames(covariates))
    assert_that(0 %in% covariates$TIME)
  } else if (length(covariates)==0) {
    covariates <- NULL
    # nothing, all good
  } else {
    stop("Covariates in wrong format")
  }
  # now covariates is NULL, or a data.frame

  # IOV: create data.frame for use as covariates
  if(length(iovParameters) > 0) {
    ## Move the IOV parameters into time-varying covariates
    if(!is.null(extraArguments$covsInterpolation) && extraArguments$covsInterpolation != "locf")
      stop("When using IOV, no interpolation for time-varying covariates should be used!")

    df <- list()
    df$OCC <- unique(regimen$OCC)
    for(i in iov) df[[i]] <- iovParameters[ names(iovParameters) == i ]

    df <- as.data.frame(df) %>% merge(regimen[, c("TIME", "OCC")])
    if(!is.null(covariates)) {
      covariates <- covariates %>%
        merge(df, by="TIME", all=T) %>% #and it is immediately sorted as well!
        dplyr::mutate_all(zoo::na.locf)
    } else {
      covariates <- df
    }
  }

  ## Make sure that all model input is defined now
  i <- modVars$params %in% c(names(parameters), colnames(covariates))
  if( any(!i) ) {
    stop("Model parameter(s) ", paste(modVars$params[i], collapse=", "), " is missing.")
  }

  # All arguments look good, let's prepare the simulation
  ev <- RxODE::eventTable()
  if(length(times) > 0) ev$add.sampling(time=times)
  ev <- addRegimenToEventTable(ev, regimen)

  # Treat missing times
  if(!is.null(covariates)) {
    missingTimes <- covariates$TIME[ !(covariates$TIME %in% ev$get.sampling()$time) ]
    if(length(missingTimes) > 0) ev$add.sampling(missingTimes)

    eventTable <- ev$get.EventTable()
    covariates <-  eventTable %>% dplyr::transmute(TIME=eventTable$time) %>% dplyr::left_join(covariates, by="TIME")
    for(i in colnames(covariates))
      covariates[, i] <- zoo::na.locf(covariates[,i]) # last observation carried forward for NA values
  }

  # Run the simulation
  result <- do.call(RxODE::rxSolve, c( list(object=model, events=ev, params=parameters, covs=covariates), extraArguments))

  # Only get the values we want
  result <- as.data.frame(result, row.names=NULL)
  names(result)[names(result)=="time"] <- "TIME"

  # Remove spurious sampling times due to covariates
  # But keep this a data.frame!
  result <- subset( result, result$TIME %in% times, drop=FALSE)

  result
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
