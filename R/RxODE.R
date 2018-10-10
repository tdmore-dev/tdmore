#' Create a TDM-capable model from an RxODE model
#'
#' @param model the RxODE model
#' @param parameters list of parameter names, or NULL to use all parameters names from RxODE
#' @param omega omega variance-covariance matrix, or NULL to use a diagonal matrix of variances 1
#' @param covariates the model covariates, named vector, or data.frame with column 'TIME', and at least TIME 0
#' @param add additive residual error, as stdev
#' @param prop proportional residual error, as stdev
#' @param exp exponential residual error, as stdev. The exponential error cannot be used in conjunction with the additive or proportional error
#' @param ... extra arguments will be passed to the model_predict call
#'
#' @return An object of class tdmore, which can be used to estimate posthoc Bayesian parameters
#' @export
#'
#' @example inst/examples/RxODE.R
tdmore.RxODE <- function(model, parameters=NULL, omega=NULL, covariates=NULL, add=0, prop=0, exp=0, ...) {
  assert_that(class(model) %in% c("RxODE")) #currently, only RxODE is supported

  # Check that parameters + covariates together supplies the parameters needed for the model
  modVars <- model$get.modelVars()
  if(is.null(parameters)) parameters=modVars$params
  if(is.null(covariates)) covariates = modVars$params[ ! modVars$params %in% parameters ] #the rest

  assert_that( all.equal( sort(c(parameters, covariates)), sort(modVars$params)) )

  if(is.null(omega)) omega = diag(rep(1, length(parameters)))
  assertthat::are_equal(nrow(omega), length(parameters))
  assertthat::are_equal(ncol(omega), length(parameters))
  assertthat::assert_that(all( eigen(omega)$values >= 0)) #matrix should be positive semi-definite

  assertthat::is.number(add)
  assertthat::is.number(prop)
  assertthat::is.number(exp)
  ## Exponential and add/prop are mutually exclusive
  if(exp != 0) assert_that(add==0 & prop==0)
  if(add != 0 || prop != 0) assert_that(exp == 0)

  structure(list(
    model=model,
    omega=omega,
    res_var=list(add=add, prop=prop, exp=exp),
    parameters=parameters,
    extraArguments=list(...)
  ), class="tdmore")
}

#' Predict for RxODE
#'
#' @param model The model itself.
#' @param newdata
#' A dataframe with at least a column 'TIME' and other values that should be observed,
#' or a numeric TIME vector to produce all values that can be observed by this model
#' @param regimen dataframe with column 'TIME' and adhering to standard NONMEM specifications otherwise (columns AMT, RATE, CMT).
#' @param parameters either a dataframe with column 'TIME' and a column for each covariate and parameter, or a named numeric vector
#' @param covariates named vector, or data.frame with column 'TIME', and at least TIME 0
#' @param extraArguments extra arguments to use
#'
#' @export
#'
#' @return
#' A data.frame similar to the observed data frame, but with predicted values.
#' @importFrom RxODE eventTable
#' @importFrom RxODE rxSolve
#' @importFrom dplyr transmute
#' @importFrom dplyr left_join
#'
model_predict.RxODE <- function(model, newdata, regimen=data.frame(TIME=c()), parameters=c(), covariates=NULL, extraArguments=list()) {
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

  # Verify arguments are good
  samplingTimes <- NULL
  oNames <- NULL
  if("data.frame" %in% class(newdata)) {
    assert_that("data.frame" %in% class(newdata))
    assert_that("TIME" %in% colnames(newdata))
    assert_that(!is.unsorted(newdata$TIME))
    oNames <- colnames(newdata)
    oNames <- oNames[oNames != "TIME"]
    assert_that(length(oNames) > 0) ## TODO: What should happen if newdata contains a column that is not predicted by the model??
    oNames <- oNames[oNames %in% c(modVars$lhs, modVars$state)] #only keep the ones that are required
    assert_that(all(oNames %in% c(modVars$lhs, modVars$state))) #TODO: not required when providing covariates! They are in rhs!
    samplingTimes <- newdata$TIME
  } else {
    assert_that(is.numeric(newdata))
    assert_that(!is.unsorted(newdata))
    samplingTimes <- newdata
    newdata <- data.frame(TIME=newdata)
  }

  assert_that("data.frame" %in% class(regimen))
  assert_that(all(c("TIME", "AMT") %in% colnames(regimen)))
  assert_that(all(colnames(regimen) %in% c("TIME", "AMT", "RATE", "DURATION", "CMT", "II", "ADDL")))
  if("II" %in% colnames(regimen) || "ADDL" %in% colnames(regimen))
    assert_that(all(c("II", "ADDL") %in% colnames(regimen)))


  # All arguments look good, let's prepare the simulation
  ev <- RxODE::eventTable()
  ev$add.sampling(time=samplingTimes)
  for(i in 1:nrow(regimen)) {
    row <- regimen[i, ,drop=FALSE]
    dosing.to = 1
    rate = NULL
    nbr.doses = 1
    dosing.interval = 24
    if(isTRUE(row$II > 0)) {
      dosing.interval <- row$II
      nbr.doses <- row$ADDL
      #duration <- max(observed$TIME) - row$TIME
      #nbr.doses <- ceiling(duration / dosing.interval)
    }

    if(isTRUE(is.finite(row$RATE))) rate=row$RATE
    if(isTRUE(is.finite(row$DURATION))) {
      if(!is.null(rate)) stop("Cannot specify RATE and DURATION in the same treatment row ", i)
      rate = row$AMT / row$DURATION
    }

    if("CMT" %in% names(row)) dosing.to <- row$CMT
    ev$add.dosing(start.time = row$TIME,
                  dose=row$AMT,
                  dosing.to=dosing.to,
                  rate=rate,
                  dosing.interval=dosing.interval,
                  nbr.doses=nbr.doses)
  }


  ## Set up the parameters
  params = NULL
  assert_that(is.numeric(parameters))
  pNames <- names(parameters)
  params = parameters

  if("data.frame" %in% class(covariates)) {
    assert_that("TIME" %in% colnames(covariates))
    assert_that(0 %in% covariates$TIME)

    missingTimes <- covariates$TIME[ !(covariates$TIME %in% ev$get.sampling()$time) ]
    if(length(missingTimes) > 0) ev$add.sampling(missingTimes)

    covs <- ev$get.EventTable() %>% dplyr::transmute(TIME=time) %>% dplyr::left_join(covariates, by="TIME")
    for(i in colnames(covs))
      covs[, i] <- zoo::na.locf( covs[, i] )

    cNames <- colnames(covariates)
    cNames <- cNames[cNames != "TIME"]
    pNames <- c(pNames, cNames)
  } else if (is.numeric(covariates)){
    pNames <- c(pNames, names(covariates))
    params = c(parameters, covariates)
    covs = NULL
  } else if (length(covariates) == 0 || is.null(covariates)) {
    covs=NULL
  } else {
    stop("Covariates in wrong format")
  }

  assert_that(all(pNames %in% modVars$params))
  assert_that(all(modVars$params %in% pNames))

  # Run the simulation
  result <- do.call(RxODE::rxSolve, c( list(object=model, events=ev, params=params, covs=covs), extraArguments))
  #result <- model$solve(events=ev, params=params, covs=covs)

  # Only get the values we want
  result <- as.data.frame(result)
  result$TIME <- result$time
  if(!is.null(oNames)) result <- result[, c("TIME", oNames)]
  result <- subset( result, TIME %in% newdata$TIME ) # Remove spurious sampling times due to covariates
  result
}
