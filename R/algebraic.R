#' Use an external function to generate an algebraic model, capable of being used with tdmore()
#' Algebraic functions are assumed to be dose-linear, meaning that they can be calculated for a single dose
#' and then summed up for multiple doses.
#'
#' @param fun a function with 'time' as the first argument, and additional parameters as needed.
#' In case parameters will be estimated, they should have a mean of '0' and be normal-distributed.
#'
#' The function can also include parameters `TIME`, `AMT`, `II`, `ADDL`, `RATE`, `DURATION`, `CMT`.
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
  pNames <- argNames[! argNames %in% c(regimenNames)]

  structure(list(
    predictFunction = function(times, regimen, parameters, covariates, iov) {
      if(!all(colnames(regimen) %in% c(regimenNames, "OCC")) || !all(regimenNames %in% colnames(regimen)))
        stop("Algebraic function requires regimen names `", paste(regimenNames, collapse=", "), "',
             but regimen provided the following columns: `", paste(colnames(regimen), collapse=", "), "'")
      if(!all(names(parameters) %in% pNames) || !all(pNames %in% c(names(parameters), names(covariates))))
        stop("Algebraic function requires parameters `", paste(pNames, collapse=", "), "',
             but was provided the following: `", paste(names(parameters), collapse=", "), "'")

      if(anyNA(regimen)) stop("The provided regimen contains NA. Cannot calculate algebraic model...")

      CONCs = apply(regimen, 1, function(regimenRow) {
        iovPrediction <- !is.null(iov)
        iovIndexes <- which(names(parameters) %in% iov)
        iovParameters <- parameters[iovIndexes]
        parameters <- parameters[!duplicated(names(parameters))]

        # Use IOV value corresponding to regimen occasion
        if(iovPrediction) {
          occasion <- regimenRow[['OCC']]
          for(iov_term in iov) {
            iovValueIndexes <- which(names(iovParameters)==iov_term)
            iovValue <- iovParameters[[iovValueIndexes[occasion]]]
            if(is.null(iovValue)) {
              stop(paste("Missing IOV values for", iov_term))
            }
            parameters[iov_term] <- iovValue
          }
        }

        # Fixed covariates
        if(is.null(covariates) || is.numeric(covariates)) {
          parameters <- c(parameters, covariates)

        # Time-varying covariates (only works at regimen times)
        } else {
          covariates <- covariates[order(covariates$TIME),]
          cNames <- colnames(covariates)
          cNames <- cNames[cNames != "TIME"]
          selectedRow <- tail(which((covariates$TIME <= as.numeric(regimenRow[['TIME']]))==TRUE), n=1)
          covs <- unlist(covariates[selectedRow, cNames])
          parameters <- parameters[!(names(parameters) %in% names(covariates))]
          parameters <- c(parameters, covs)
        }

        args <- list(times)
        names(args) <- tArg

        args <- c(args, as.list(regimenRow[!(names(regimenRow) %in% c("OCC"))])) # add regimen names
        args <- c(args, parameters) # add parameters

        funValue <- do.call(fun, args=args)
        ifelse(times < args$TIME, 0, funValue)
      })
      if(length(times)==0) return(numeric())
      if(length(times)==1) return( sum(CONCs, na.rm=TRUE) )
      if(nrow(regimen) == 0) return( rep(0, length.out=length(times)) )

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
#' @param covariates named numeric vector, or data.frame with column 'TIME', and at least TIME 0
#' @param iov character array with the IOV terms, NULL if no IOV, IOV on KA only
#'
#' @return
#' A data.frame similar to the newdata data frame, but with CONC column filled out.
#'
#' @export
#' @keywords internal
model_predict.algebraic <- function(model, times, regimen=data.frame(TIME=numeric(), AMT=numeric(), OCC=numeric()), parameters=numeric(), covariates=NULL, iov=NULL, extraArguments=list()) {

  # Check times and regimen objects
  checkTimes(times)
  checkRegimen(regimen, iov)
  regimen <- flatten(regimen)

  # Check parameters and covariates
  assertthat::assert_that(is.numeric(parameters))

  # Check covariates
  if (is.null(covariates) || is.numeric(covariates)) {
    pNames <- names(c(parameters, covariates))
  } else if (is.data.frame(covariates)) {
    cNames <- colnames(covariates)
    cNames <- cNames[cNames != "TIME"]
    pNames <- unique(c(names(parameters), cNames))
  } else {
    stop("Covariates should be either a numeric vector or data frame")
  }
  assertthat::assert_that(all(pNames %in% model$parameters))
  assertthat::assert_that(all(model$parameters %in% pNames))

  # Predict concentrations
  conc <- model$predictFunction(times, regimen, parameters, covariates, iov)

  return(data.frame(TIME=times, CONC=conc))
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
