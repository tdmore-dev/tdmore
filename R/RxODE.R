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


#' Prepare a cache object for RxODE simulation
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
#' A list, to be used as a cache object for model_predict
#'
#' @importFrom RxODE eventTable rxSolve
#' @importFrom dplyr transmute left_join bind_rows arrange mutate_all
#'
#' @keywords internal
#'
model_prepare.RxODE <- function(model, times, regimen=data.frame(TIME=numeric()), parameters=numeric(), covariates=NULL, iov=NULL, extraArguments=list()) {
  modVars <- model$get.modelVars()
  # RxODE does not allow to simulate 'nothing'
  # Manually construct an empty data.frame with the right columns
  if(length(times) == 0) {
    colNames <- c("TIME", modVars$lhs, modVars$state)
    df <- data.frame(matrix(numeric(), nrow=0, ncol=length(colNames),
                            dimnames=list(c(), colNames)),
                     stringsAsFactors=F)
    cache <- list(output=df)
    return(cache)
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
    stop("Model parameter(s) ", paste(modVars$params[!i], collapse=", "), " is missing.")
  }

  # All arguments look good, let's prepare the simulation
  ev <- RxODE::eventTable()
  if(length(times) > 0) ev$add.sampling(time=times)
  ev <- addRegimenToEventTable(ev, regimen, extraArguments$nbSSDoses)

  # Treat missing times
  if(!is.null(covariates)) {
    missingTimes <- covariates$TIME[ !(covariates$TIME %in% ev$get.sampling()$time) ]
    if(length(missingTimes) > 0) ev$add.sampling(missingTimes)

    eventTable <- ev$get.EventTable()
    covariates <-  eventTable %>% dplyr::transmute(TIME=eventTable$time) %>% dplyr::left_join(covariates, by="TIME")
    for(i in colnames(covariates))
      covariates[, i] <- zoo::na.locf(covariates[,i]) # last observation carried forward for NA values
  }

  cache <- list()
  cache$parameters <- function(x) {
    x <- x[ ! names(x) %in% iov ]
    parameters[names(x)] <- x
    if(length(parameters) == 0) return(NULL)
    parameters
  }

  update_table <- function(dt, x) {
    ### We assume the order of parameters matches the IOV order
    x <- x[ names(x) %in% iov ]
    if(length(x)==0) return(dt)

    ## Update the IOV terms in the covs data.frame
    NOcc <- length(x) / length(iov)
    Niov <- length(iov)
    for(i in seq_len(NOcc)) {
      thisOcc <- seq_len(Niov) + (i-1)*Niov
      for(j in seq_along(iov)) {
        data.table::set(dt, which(dt$OCC==i), iov[j], x[ (i-1)*Niov + j] )
      }
    }
  }

  cutoffVersion <- base::package_version("0.8.1")
  rxVersion <- utils::packageVersion('RxODE')
  if(rxVersion >= cutoffVersion) {
    ## new version of RxODE
    if(!is.null(covariates)) ev <- cbind(ev, covariates)
    ev <- data.table::as.data.table(ev)
    cache$ev <- function(x) {
      update_table(ev, x)
      ev
    }
    cache$rxSolveArgs <- function(parameters) {
      list(
        events=cache$ev(parameters),
        params=cache$parameters(parameters)
        ## covs cannot even be specified anymore!!
      )
    }
  } else {
    ## old version, covariates via 'covs' argument
    if(!is.null(covariates)) covariates <- data.table::as.data.table(covariates)
    cache$covs <- function(x) {
      if(is.null(covariates)) return(NULL)
      update_table(covariates, x)
      covariates
    }
    cache$ev <- function(x) { ev }
    cache$rxSolveArgs <- function(parameters) {
      list(
        events=cache$ev(parameters),
        params=cache$parameters(parameters),
        covs=cache$covs(parameters)
      )
    }
  }

  return(cache)
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
#'
#' @keywords internal
#'
model_predict.RxODE <- function(model, times, regimen=data.frame(TIME=numeric()), parameters=numeric(), covariates=NULL, iov=NULL, extraArguments=list(), cache=NULL) {
  if(is.null(cache)) {
    cache <- model_prepare(model, times, regimen, parameters, covariates, iov, extraArguments)
  }

  if(!is.null( cache$output) ) return(cache$output) ## output always same, no matter the parameters

  # Run the simulation
  result <- do.call(RxODE::rxSolve, c( list(object=model), cache$rxSolveArgs(parameters), extraArguments))

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
#' @param nbSSDoses number of doses to prepend for simulation of steady-state
#'
#' @return a completed event table
#'
addRegimenToEventTable <- function(eventTable, regimen, nbSSDoses=NULL) {
  if(is.null(nbSSDoses)) nbSSDoses=5 #default
  for(i in seq_len(nrow(regimen))) {
    row <- regimen[i, ,drop=FALSE]
    dosing.to = 1
    rate = NULL
    nbr.doses = 1
    dosing.interval = row$II
    start.time = row$TIME

    if("ADDL" %in% colnames(row)) {
      if(is.null(dosing.interval)) stop("Please define the dosing interval II in order to use ADDL doses")
      nbr.doses <- nbr.doses + row$ADDL
    }
    if (isTRUE(row$SS == 1)) {
      if(is.null(dosing.interval)) stop("Please define the dosing interval II in order to use SS=1")

      ## TODO: Use new functionality of new RxODE??
      nbr.doses=nbr.doses + nbSSDoses
      start.time = start.time - row$II * nbSSDoses
      if(any( eventTable$get.dosing()$time > start.time ))
        warning("Possible collision of steady-state dose on ", row$TIME, " with other treatments...")
    }

    if("RATE" %in% names(row) && isTRUE(is.finite(row$RATE))) rate=row$RATE
    if("DURATION" %in% names(row) && isTRUE(is.finite(row$DURATION))) {
      if(!is.null(rate)) stop("Cannot specify RATE and DURATION in the same treatment row ", i)
      rate = row$AMT / row$DURATION
    }

    #strange defaults in RxODE
    #see https://github.com/nlmixrdevelopment/RxODE/blob/68ecc64b0fd7e231ac0ff38541715bbc8031f583/src/et.cpp#L2485
    if(nbr.doses == 1) dosing.interval <- 24

    if("CMT" %in% names(row)) dosing.to <- row$CMT
    eventTable$add.dosing(start.time = start.time,
                  dose=row$AMT,
                  dosing.to=dosing.to,
                  rate=rate,
                  dosing.interval=dosing.interval,
                  nbr.doses=nbr.doses)
  }
  return(eventTable)
}
