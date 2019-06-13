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
#' The function returns a single vector with the concentrations for the given times.
#'
#' It can also return a named list of values.
#' The vector with the output name will be summed up across multiple administrations.
#' Any other values will be treated as parameter values. The output will be observed using last-observation-carry-forward.
#'
#' @param output Name for the output
#'
#' @return An algebraic prediction model
#' @importFrom utils tail
#' @export
algebraic <- function(fun, output="CONC") {
  if(length(output) != 1) stop("You should only specify a single 'output'")
  argNames <- names(formals(fun))

  tArg <- argNames[1]
  argNames <- argNames[-1] #remove first argument; it is always the evaluation time

  regimenNames <- c("TIME", "AMT", "II", "ADDL", "RATE", "DURATION", "CMT", "SS")
  regimenNames <- argNames[argNames %in% regimenNames]
  pNames <- argNames[! argNames %in% c(regimenNames)]

  structure(list(
    output=output,
    predictFunction = function(times, regimen, parameters, covariates, iov) {
      if(!all(colnames(regimen) %in% c(regimenNames, "OCC")) || !all(regimenNames %in% colnames(regimen)))
        stop("Algebraic function requires regimen names `", paste(regimenNames, collapse=", "), "',
             but regimen provided the following columns: `", paste(colnames(regimen), collapse=", "), "'")

      if(anyNA(regimen)) stop("The provided regimen contains NA. Cannot calculate algebraic model...")

      if(length(times) ==0) {
        df <- data.frame(TIME=numeric(0))
        df[, output] <- numeric(0)
        return(df)
      }

      regimen <- regimen[ regimen$TIME <= max(times, -Inf), ]

      df <- data.frame(TIME=times)
      if(nrow(regimen) == 0) {
        df[, output] <- rep(0, length.out=length(times))
        return( df ) #not entirely happy, because the extra output will not be in there...
      }

      regimenResult = apply(regimen, 1, function(regimenRow) {
        iovPrediction <- !is.null(iov)
        iovIndexes <- which(names(parameters) %in% iov)
        iovParameters <- parameters[iovIndexes]
        parameters <- parameters[!duplicated(names(parameters))]

        # Use IOV value corresponding to regimen occasion
        if(iovPrediction) {
          occasion <- regimenRow[['OCC']]
          for(iov_term in iov) {
            iovValueIndexes <- which(names(iovParameters)==iov_term)
            if(occasion > length(iovValueIndexes))
              stop(paste("Missing IOV values for", iov_term))
            iovValue <- iovParameters[[iovValueIndexes[occasion]]]
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
          selectedRow <- utils::tail(which((covariates$TIME <= as.numeric(regimenRow[['TIME']]))==TRUE), n=1)
          covs <- unlist(covariates[selectedRow, cNames])
          parameters <- parameters[!(names(parameters) %in% names(covariates))]
          parameters <- c(parameters, covs)
        }

        args <- list(times)
        names(args) <- tArg

        args <- c(args, as.list(regimenRow[!(names(regimenRow) %in% c("OCC"))])) # add regimen names
        args <- c(args, parameters) # add parameters

        funValue <- do.call(fun, args=args)
        if(!is.list(funValue)) {
          funValue[times < args$TIME] <- 0
          out <- list(funValue)
          names(out) <- output
          out
        } else {
          stopifnot( output %in% names(funValue))
          funValue[[output]][times < args$TIME] <- 0

          n <- lapply(funValue, length)
          stopifnot( n[[output]] == length(times) )
          stopifnot( all( n[ names(n) != output ] == 1 ))
          funValue
        }
      })
      #We now have a list of all results from each treatment

      ## We need to SUM the concentrations
      ## And fetch the appropriate parameter values from the right regimen result
      CONCs <- vapply(regimenResult, function(x) { x[[output]] }, FUN.VALUE=numeric(length(times)) )

      if(length(times)==0) {
        for(outputCol in names( regimenResult[[1]] ))df[, outputCol] <- numeric()
      } else if (length(times)==1) {
        df[, output] <- sum(CONCs)
        i <- max(which(times >= regimen$TIME))
        otherNames <- setdiff( names(regimenResult[[1]]), output)
        for(j in otherNames) df[, j] <- vapply( regimenResult[i], function(x){ x[[j]] }, numeric(1))
      } else {
        df[, output] <- apply(CONCs, 1, sum) ## why was na.rm here?
        i <- vapply(times, function(t) {
          max(which(t >= regimen$TIME))
        }, integer(1)) #find the nearest treatment regimen
        otherNames <- setdiff( names(regimenResult[[1]]), output)
        for(j in otherNames) df[, j] <- vapply( regimenResult[i], function(x){ x[[j]] }, numeric(1))
      }
      return(df)
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
model_predict.algebraic <- function(model, times, regimen=NULL, parameters=numeric(), covariates=NULL, iov=NULL, extraArguments=list(), cache=NULL) {
  if(is.null(regimen)){
    regimen <- data.frame(TIME=numeric(), AMT=numeric())
    if(!is.null(iov)) regimen$OCC <- numeric()
  }

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

  if(!setequal(pNames, model$parameters))
    stop("Algebraic function requires parameters `", paste(model$parameters, collapse=", "), "',
             but was provided the following: `", paste(pNames, collapse=", "), "'")

  # Predict concentrations
  model$predictFunction(times, regimen, parameters, covariates, iov)
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
