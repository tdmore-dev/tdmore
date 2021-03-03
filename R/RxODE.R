#' @name tdmore
#'
#' @param res_var the residual variability
#' @param parameters list of parameter names, or NULL to automatically detect
#' The automatic detection will analyze omega first, to see if there are names present.
#' If not, it will use all parameters from the RxODE model
#' @param omega omega variance-covariance matrix, or NULL to use a diagonal matrix of variance 1 for all input parameters
#' @param iov list of parameter names related to IOV, NULL if no IOV
#' @param ... extra arguments will be passed to the model_predict call
#'
#' @note You can use a named omega parameter to distinguish between unexplained variability (described by an a priori distribution) and
#' a covariate (assumed to be known for all individuals). Any input parameters not provided in omega are assumed to be covariates.
#'
#' @export
#'
#' @example inst/examples/RxODE.R
tdmore.RxODE <- function(model, res_var, parameters=NULL, omega=NULL, iov=NULL, ...) {
  if(is.null(parameters)) {
    # try to guess using the omega matrix
    if(!is.null(names(omega)) ) parameters <- names(omega)
    else if(!is.null(rownames(omega)) ) parameters <- rownames(omega)
    else if(!is.null(colnames(omega)) ) parameters <- colnames(omega)
    else {
      # keep it NULL; get the parameters from RxODE model instead
    }
  }
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
#' @noRd
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
  # XXX: This is not necessary anymore
  #if(!RxODE::rxDllLoaded(model)) {
  #  warning("RxDLL was not loaded, loading...")
  #  RxODE::rxLoad(model)
  #  if(!RxODE::rxDllLoaded(model)) stop("Tried to reload RxDLL, but failed. Cannot continue!")
  #}
  if(!is.null(iov) && length(iov) > 0 && !"OCC" %in% colnames(regimen)) {
    shouldWarn <- getOption("tdmore.warnIov", NA)
    if(is.na(shouldWarn)) {
      warning("Adding OCC column to regimen\nThis warning will appear once during the runtime of this application. To enable consistently, set options(tdmore.warnIov)")
      options(tdmore.warnIov=FALSE)
    } else if (shouldWarn) {
      warning("Adding OCC column to regimen")
    } else {
      # do not warn
    }
    regimen$OCC <- seq_len(nrow(regimen))
  }

  # Check times and regimen objects
  checkTimes(times)
  checkRegimen(regimen, iov)

  # Check parameters
  stopifnot(is.numeric(parameters))
  i <- names(parameters) %in% iov
  iovParameters <- parameters[i] #split the parameters in two arrays
  parameters <- parameters[!i]

  # Set up the covariates
  if(is.data.frame(covariates) && nrow(covariates)==0) {
    #ignore; treat it as NULL
    covariates <- NULL
  } else if (is.numeric(covariates)){
    i <- is.na(covariates)
    if(any(i)) stop("Covariate ", names(covariates)[i], " has NA value, cannot predict")
    parameters <- c(parameters, covariates)
    covariates <- NULL
  } else if(is.data.frame(covariates)) {
    stopifnot("TIME" %in% colnames(covariates))
    stopifnot(0 %in% covariates$TIME)
  } else if (length(covariates)==0) {
    covariates <- NULL
    # nothing, all good
  } else {
    stop("Covariates in wrong format")
  }
  # now covariates is NULL, or a data.frame

  # All arguments look good, let's prepare the simulation
  regimen$AMT[is.na(regimen$AMT)] <- 0 #remove NA
  ev <- data.table::rbindlist( c(list(list(time=times, evid=0)), as.RxODE_regimen(regimen)), fill=TRUE ) %>%
    dplyr::arrange(.data$time)

  # Add the covariates into EV
  if(!is.null(covariates)) {
    covariates$evid <- 2
    covariates$time <- covariates$TIME
    covariates <- covariates[, setdiff(names(covariates), "TIME") ]
    ev <- data.table::rbindlist(list(covariates, ev), fill=TRUE) %>% dplyr::arrange(.data$time)
    for(i in setdiff(names(covariates), c("time", "evid"))) {
      value <- zoo::na.locf(ev[, i], na.rm=FALSE)
      if(is.na(value[1])) stop("Covariate ", i, " has NA starting value, cannot predict")
      ev[,i] <- value
    }
  }

  ev <- ev %>% dplyr::arrange(.data$time, factor(.data$evid, levels=c(1, 0, 2) )) #FIRST regimen, so OCC can propagate correctly
  if("OCC" %in% names(ev)) ev$OCC <- zoo::na.locf(ev$OCC, na.rm=FALSE)
  if("OCC" %in% names(ev)) ev$OCC[is.na(ev$OCC)] <- 1 #leading NA should be '1'
  ev <- ev %>% dplyr::arrange(.data$time, factor(.data$evid, levels=c(2, 0, 1) )) #FIRST covariate, then OBS then TMT

  if(length(iovParameters) > 0) {
    ## Move the IOV parameters into the EV data.frame
    if(!is.null(extraArguments$covsInterpolation) && extraArguments$covsInterpolation != "locf")
      stop("When using IOV, no interpolation for time-varying covariates should be used!")

    df <- list()
    df$OCC <- unique(regimen$OCC)
    N <- length(df$OCC)
    for(i in iov) {
      par <- iovParameters[ names(iovParameters) == i ]
      par <- c(par, rep(0, N-length(par))) #pad with extra zero for new occasions
      df[[i]] <- par
    }

    # Add parameters to ev
    ev <- dplyr::left_join(ev, as.data.frame(df), by="OCC")
  }

  ## Make sure that all model input is defined now
  i <- modVars$params %in% c(names(parameters), colnames(ev))
  if( any(!i) ) {
    stop("Model parameter(s) ", paste(modVars$params[!i], collapse=", "), " is missing.")
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
#' @noRd
model_predict.RxODE <- function(model, times, regimen=data.frame(TIME=numeric()), parameters=numeric(), covariates=NULL, iov=NULL, extraArguments=list(), cache=NULL) {
  if(is.null(cache)) {
    cache <- model_prepare(model, times, regimen, parameters, covariates, iov, extraArguments)
  }

  if(!is.null( cache$output) ) return(cache$output) ## output always same, no matter the parameters

  # Run the simulation
  result <- do.call(RxODE::rxSolve, c( list(object=model, returnType="data.frame", addDosing=NULL), cache$rxSolveArgs(parameters), extraArguments))
  names(result)[names(result)=="time"] <- "TIME"

  # Remove spurious sampling times due to covariates
  # But keep this a data.frame (drop=FALSE)
  result <- subset( result, result$TIME %in% times, drop=FALSE)

  result
}

#' Add the given regimen into the RxODE event table.
#'
#' @param regimen the specified regimen
#' @param nbSSDoses number of doses to prepend for simulation of steady-state
#'
#' @return a completed event table
#' @noRd
as.RxODE_regimen <- function(regimen, nbSSDoses=NULL) {
  result <- list()
  for(i in seq_len(nrow(regimen))) {
    row <- regimen[i, ,drop=FALSE]
    args <- list()
    args$time <- row$TIME
    args$amt = row$AMT
    if("II" %in% names(row)) args$ii <- row$II
    if("CMT" %in% names(row)) args$cmt <- row$CMT

    if("ADDL" %in% colnames(row)) {
      if(is.null(args$ii)) stop("Please define the dosing interval II in order to use ADDL doses")
      args$addl <- row$ADDL
    }
    if (isTRUE("SS" %in% names(row) && row$SS > 0)) {
      if(is.null(args$ii)) stop("Please define the dosing interval II in order to use SS=1")

      if(is.null(nbSSDoses)) {
        args$ss <- row$SS
      } else {
        if(is.null(args$addl)) args$addl <- 0
        args$addl = args$addl + nbSSDoses
        args$time = args$time - row$II * nbSSDoses
        # if(any( eventTable$get.dosing()$time > args$time ))
        #   warning("Possible collision of steady-state dose on ", row$TIME, " with other treatments...")
      }
    }

    if("ii" %in% names(args) && ! any(c("ss", "addl") %in% names(args) ) ) {
      #RxODE is not happy with II and no additional doses or steady state dosing
      #let's simply remove II
      args$ii <- NULL
    }

    if("RATE" %in% names(row) && isTRUE(is.finite(row$RATE))) args$rate=row$RATE
    if("DURATION" %in% names(row) && isTRUE(is.finite(row$DURATION))) {
      if(!is.null(args$rate)) stop("Cannot specify RATE and DURATION in the same treatment row ", i)
      args$dur = row$DURATION
    }

    args$evid=1 #dose
    if("OCC" %in% names(row)) args$OCC = row$OCC
    result[[i]] <- args
  }
  return(result)
}
