#' Use an external function to generate an algebraic model, capable of being used with tdmore()
#' Algebraic functions are assumed to be dose-linear, meaning that they can be calculated for a single dose
#' and then summed up for multiple doses.
#'
#' @param fun a function with 'time' as the first argument, and additional parameters as needed.
#' In case parameters will be estimated, they should have a mean of '0' and be normal-distributed.
#'
#' The function can also include parameters `TIME`, `AMT`, `II`, `ADDL`, `RATE`, `DURATION` or `CMT`.
#' These are then taken from the treatment regimen.
#'
#' @return An algebraic prediction model
#' @export
algebraic <- function(fun) {
  argNames <- names(formals(fun))

  tArg <- argNames[1]
  argNames <- argNames[-1] #remove first argument; it is always the evaluation time

  regimenNames <- c("TIME", "AMT", "II", "ADDL", "RATE", "DURATION", "CMT")
  regimenNames <- argNames[argNames %in% regimenNames]
  pNames <- argNames[! argNames %in% regimenNames]

  structure(list(
    predictFunction = function(times, regimen, parameters) {
      if(!all(colnames(regimen) %in% regimenNames) || !all(regimenNames %in% colnames(regimen)))
        stop("Algebraic function requires regimen names `", paste(regimenNames, collapse=", "), "',
             but regimen provided the following columns: `", paste(colnames(regimen), collapse=", "), "'")
      if(!all(names(parameters) %in% pNames) || !all(pNames %in% names(parameters)))
        stop("Algebraic function requires parameters `", paste(pNames, collapse=", "), "',
             but was provided the following: `", paste(names(parameters), collapse=", "), "'")
      CONCs = apply(regimen, 1, function(x) {
        args <- list(times)
        names(args) <- tArg

        args <- c(args, as.list(x)) # add regimen names

        args <- c(args, parameters) # add parameters
        do.call(fun, args=args)
      })
      CONC <- apply(CONCs, 1, sum, na.rm=TRUE)
      return(CONC)
    },
    parameters = pNames
  ),
  class = "algebraic")
}

#' Predict for algebraic models
#'
#' We only support dosing into the default compartment, and only bolus doses
#'
#' @param model the algebraic model
#' @param newdata dataframe with a column 'TIME' and 'CONC'. Other columns are ignored.
#' @param regimen dataframe with column 'TIME' and adhering to standard NONMEM specifications otherwise (columns AMT, RATE, CMT)
#' @param parameters a named numeric vector
#' @param covariates named vector, or data.frame with column 'TIME', and at least TIME 0
#'
#' @return
#' A data.frame similar to the newdata data frame, but with CONC column filled out.
#'
#' @export
#' @keywords internal
model_predict.algebraic <- function(model, newdata, regimen=data.frame(TIME=c()), parameters=c(),
                                    covariates=NULL, extraArguments=list()) {
  # Verify arguments are good
  assertthat::assert_that("data.frame" %in% class(newdata))
  assertthat::assert_that("TIME" %in% colnames(newdata))
  oNames <- colnames(newdata)
  oNames <- oNames[oNames != "TIME"]
  assertthat::assert_that(length(oNames) > 0)
  assertthat::assert_that(all(oNames %in% c("CONC")))

  assertthat::assert_that("data.frame" %in% class(regimen))
  assertthat::assert_that(all(c("TIME", "AMT") %in% colnames(regimen)))
  assertthat::assert_that(all(colnames(regimen) %in% c("TIME", "AMT", "RATE", "II", "ADDL")))
  # if either II or ADDL is mentioned, the other one needs to be present as well
  if("II" %in% colnames(regimen) || "ADDL" %in% colnames(regimen))
    assertthat::assert_that(all(c("II", "ADDL") %in% colnames(regimen)))

  assertthat::assert_that(is.numeric(parameters))
  pNames <- names(parameters)
  params = parameters
  assertthat::assert_that(all(pNames %in% model$parameters))
  assertthat::assert_that(all(model$parameters %in% pNames))

  # All arguments look good, let's prepare the simulation
  if(is.data.frame(newdata)) times <- newdata$TIME
  assertthat::assert_that(is.numeric(times))

  if(is.data.frame(covariates)) stop("Time-varying covariates not supported")
  params <- c(params, covariates)

  CONC <- model$predictFunction(times, regimen, params)

  # Only get the values we want
  result <- data.frame(TIME=times, CONC=CONC)
  result[, colnames(newdata)]
}


#' Build a tdmore object based on an algebraic model
#'
#' @param model The algebraic model
#' @param res_var the residual variability
#' @param omega the omega matrix of the model
#' @param ... extra arguments will be passed to the model_predict call
#'
#' @return a tdmore object, capable of estimating bayesian individual parameters
#' @export
tdmore.algebraic <- function(model, res_var, omega, ...) {
  if(is.numeric(omega)) omega <- vectorToDiagonalMatrix(omega)
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
    covariates=covariates
  ), class="tdmore")

  # Check consistency and return
  return(checkTdmore(tdmore))
}
