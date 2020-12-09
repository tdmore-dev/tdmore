#' Initialize a structural model using algebraic equations.
#'
#' @note
#' Algebraic functions are assumed to be dose-linear, meaning that they can be calculated for a single dose
#' and then summed up for multiple doses.
#'
#' @param fun function with named arguments. The first argument is always `time` (a numeric vector), the other arguments
#' may be the initial state, the regimen (TIME, AMT, II, ADDL, RATE, DURATION, CMT, SS), or
#' the covariates or parameters (remaining arguments). These are all single numbers.
#'
#' The initial values should be viewed as 'initial value at dosing time `TIME`.
#' `AMT` may be NA, which is seen as a change in covariates (rather than a new administration).
#'
#' The function returns a named list with numeric vectors of size `time` or size `1`.
#' @param inits named numeric vector with the initial states for the model
#' These states should also be present in the output list of the function.
#' @param na.rm time-varying covariates result in a new call to the prediction function,
#' using AMT=NA. If TRUE, any AMT=NA is changed by AMT=0. Some models apply special rules
#' when a new dose is applied. They may set this switch to FALSE to receive AMT=NA for when a time-varying covariate changes
#' (but no actual treatment occurs).
#' @return An algebraic prediction model
#' @export
algebraic <- function(fun, inits=NULL, na.rm=TRUE) {
  stopifnot(is.logical(na.rm))

  argNames <- names(formals(fun))

  tArg <- argNames[1]
  argNames <- argNames[-1] #remove first argument; it is always the evaluation time

  regimenNames <- c("TIME", "AMT", "II", "ADDL", "RATE", "DURATION", "CMT", "SS")
  regimenNames <- argNames[argNames %in% regimenNames]
  pNames <- argNames[! argNames %in% c(regimenNames)]
  initsAuto <- FALSE
  if(is.null(inits)) {
    #auto-detect inits
    initsAuto <- TRUE
    definedArgs <- pNames[ ! sapply( formals(fun)[pNames] , is.name) ]
    inits <- unlist( formals(fun)[definedArgs] )
  }
  if( ! all(names(inits) %in% pNames) ) stop("Not all initial values are function arguments to `fun`")
  pNames <- setdiff(pNames, names(inits)) #remove inits

  structure(list(
    output=output,
    fun=fun,
    inits=inits,
    initsAuto=initsAuto,
    na.rm=na.rm,
    parameters=pNames
  ), class = "algebraic")
}

#' Predict for algebraic models
#'
#' @param model the algebraic model
#' @param times numeric vector at which to provide prediction values
#' @param regimen dataframe with column 'TIME' and adhering to standard NONMEM specifications otherwise (columns AMT, RATE, CMT)
#' @param parameters a named numeric vector
#' @param covariates named numeric vector, or data.frame with column 'TIME', and at least TIME 0
#' @param iov character array with the IOV terms, NULL if no IOV. These values should then be present in the parameters vector, and be repeated for
#' max(regimen$OCC).
#'
#' @return
#' A data.frame containing all output of the prediction function of the algebraic model
#'
#' @engine
#' @keywords internal
model_predict.algebraic <- function(model, times, regimen=NULL, parameters=numeric(), covariates=NULL, iov=NULL, extraArguments=list(), cache=NULL) {
  if(is.null(regimen)){
    regimen <- data.frame(TIME=numeric(), AMT=numeric())
    if(!is.null(iov)) regimen$OCC <- numeric()
  }
  # covariates is either NULL, a numeric vector, or a data.frame

  ## Construct a single 'event' data.frame and
  ## slice into time windows
  if(! 0 %in% regimen$TIME ) {
    startRegimen <- data.frame(TIME=0)
    if(!is.null(iov)) startRegimen$OCC = 1
    regimen <- dplyr::bind_rows(startRegimen, regimen)
  }

  if(is.data.frame(covariates)) {
    events <- dplyr::full_join(regimen, covariates, by="TIME")
    events <- dplyr::arrange(events, .data$TIME)
    if("OCC" %in% colnames(regimen)) events$OCC <- zoo::na.locf(events$OCC) #cover the holes
  } else {
    if(is.numeric(covariates) && length(covariates) >= 1) {
      events <- cbind(regimen, as.list(covariates))
    } else if(is.null(covariates)) {
      events <- regimen
    } else {
      stop("Invalid covariates provided: ", covariates)
    }
  }

  if(is.null(iov)) {
    if(length(parameters) > 0) events <- cbind(events, as.list(parameters))
  } else {
    ## ensure we start with occ=1
    stopifnot(events$OCC[1] == 1)
    ## ensure there are enough iovParameters for the number of occasions
    maxOcc <- max(events$OCC)
    iovParameters <- list()
    for(i in iov) iovParameters[[i]] <- parameters[ names(parameters) == i ]
    lapply(iovParameters, function(x) {stopifnot(length(x) == maxOcc) })

    # add to events list
    iovParameters <- tibble::as_tibble(iovParameters)
    iovParameters$OCC <- seq_len(maxOcc)
    events <- dplyr::left_join( events, iovParameters, by="OCC")
    events <- cbind(events, as.list(parameters[setdiff(names(parameters), iov) ]))
  }

  if(model$na.rm) events$AMT[ is.na(events$AMT) ] <- 0
  fakeDuration <- c( events$TIME[-nrow(events)] - events$TIME[-1], 1 )
  if(model$na.rm & "DURATION" %in% colnames(events)) events$DURATION[is.na(events$DURATION)] <- fakeDuration[is.na(events$DURATION)]
  if(model$na.rm & "RATE" %in% colnames(events)) events$RATE[is.na(events$RATE)] <- 0

  inits <- model$inits
  result <- list()
  for(i in seq_along(events$TIME)) {
    tStart <- events$TIME[i]
    tEnd <- if(i==nrow(events)) Inf else events$TIME[i+1]
    tObs <- times[ times >= tStart & times < tEnd ]

    funtimes <- c(tObs, tEnd)
    thisEvents <- as.list(
      events[i,
             intersect(colnames(events), names(formals(model$fun))),
             drop=F]) #only take the event names that actually feature as function arguments (ignore FIX, OCC, ...)
    arglist <- c(
      list(funtimes),
      thisEvents,
      as.list(inits)
    )

    out <- do.call(model$fun, arglist)
    out <- cbind(TIME=funtimes, as.data.frame(out))

    result[[i]] <- out[-nrow(out), ,drop=FALSE]
    availableInits <- intersect(names(inits), colnames(out)) #we only use the available initial estimates
    inits <- out[nrow(out), availableInits, drop=TRUE]
    names(inits) <- availableInits
  }

  dplyr::bind_rows(result)
}

#' Build a tdmore object based on an algebraic model
#'
#' @param model The algebraic model
#' @param res_var the residual variability
#' @param omega the omega matrix of the model
#' @param iov iov terms
#' @param ... extra arguments will be passed to the model_predict call
#'
#' @return a tdmore object, capable of estimating bayesian individual parameters
#' @export
tdmore.algebraic <- function(model, res_var, omega, iov=NULL, ...) {
  if(is.numeric(omega) && !is.matrix(omega)) omega <- vectorToDiagonalMatrix(omega)
  if(!all(colnames(omega) %in% model$parameters))
    stop("Not all omega parameters are known to the model...")

  parameters <- colnames(omega) #the OMEGA values are the parameters
  covariates <- model$parameters[ ! model$parameters %in% parameters ] #the rest are covariates

  if(class(res_var) != "list") res_var <- list(res_var)

  tdmore <- structure(list(
    model=model,
    omega=omega,
    res_var=res_var,
    parameters=parameters,
    covariates=covariates,
    iov=iov
  ), class="tdmore")

  # Check consistency and return
  return(checkTdmore(tdmore))
}
