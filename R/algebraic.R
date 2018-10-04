## Library of functions
#' A function to predict concentration in a 1cpt PK model with bolus administration
#'
#' @param THETA a list(V=xx, CL=yy) with all typical values
#' @param OMEGA a list(V=0.20, CL=0.30) with all standard deviations for IIV
#'
#' @return A function to predict PK 1cpt concentration, useable in algebraic()
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' predictFunction <- pk1cptivbolusVCL(THETA=list(V=10, CL=5), OMEGA=list(V=0.20, CL=0.30))
#' model <- algebraic(predictFunction)
#' tdmore <- tdmore(model)
pk1cptivbolusVCL <- function(THETA=list(V=10, CL=5), OMEGA=list(V=0.20, CL=0.30)) {
  function(times, regimen, ETA_V, ETA_CL) {
    assertthat::assert_that(all(colnames(regimen) == c("TIME", "AMT")))
    V = THETA$V * exp(OMEGA$V * ETA_V)
    CL = THETA$CL * exp(OMEGA$CL * ETA_CL)
    k = CL / V
    t = times

    CONC <- rep(0, length(times))
    for(i in seq(1, nrow(regimen))) {
      tD = regimen$TIME[i]
      D = regimen$AMT[i]
      II = regimen$II[i] # a number, 0, or NULL
      if(II > 0) { # Keep repeating the dose until we are past the last observation time
        nbrDoses <- row$ADDL
        while(tD <= max(t) && nbrDoses > 0) {
          CONC <- CONC + ifelse(t >= tD, D / V * exp(-k*(t-tD)), 0)
          tD <- tD + II
          nbrDoses <- nbrDoses - 1
        }
      } else {
        # Single administration
        CONC <- CONC + ifelse(t >= tD, D / V * exp(-k*(t-tD)), 0)
      }
    }
    CONC
  }
}


#' A function to predict concentration in a 1cpt PK model with bolus administration
#'
#' @param THETA a list(KA=xx, V=xx, CL=yy) with all typical values
#' @param OMEGA a list(KA=xx, V=0.20, CL=0.30) with all standard deviations for IIV
#'
#' @return A function to predict PK 1cpt concentration, useable in algebraic()
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' predictFunction <- pk1cptoralbolusVCL(THETA=list(KA=0.5, V=10, CL=5), OMEGA=list(KA=0.70, V=0.20, CL=0.30))
#' model <- algebraic(predictFunction)
#' tdmore <- tdmore(model)
pk1cptoralbolusVCL <- function(THETA=list(KA=0.5, V=10, CL=5), OMEGA=list(KA=0, V=0.20, CL=0.30)) {
  function(times, regimen, ETA_V, ETA_CL, ETA_KA) {
    assertthat::assert_that(all(colnames(regimen) == c("TIME", "AMT")))
    V = THETA$V * exp(OMEGA$V * ETA_V)
    CL = THETA$CL * exp(OMEGA$CL * ETA_CL)
    k = CL / V
    KA = THETA$KA * exp(OMEGA$KA * ETA_KA)
    t = times

    CONC <- rep(0, length(times))
    for(i in seq(1, nrow(regimen))) {
      tD = regimen$TIME[i]
      D = regimen$AMT[i]
      II = regimen$II[i] # a number, 0, or NULL
      if(II > 0) { # Keep repeating the dose until we are past the last observation time
        nbrDoses <- row$ADDL
        while(tD <= max(t) && nbrDoses > 0) {
          CONC <- CONC + ifelse(t >= tD, D / V * exp(-k*(t-tD)), 0)
          tD <- tD + II
          nbrDoses <- nbrDoses - 1
        }
      } else {
        # Single administration
        CONC <- CONC + ifelse(t >= tD, D / V * exp(-k*(t-tD)), 0)
      }
    }
    CONC
  }
}


#' Use an external function to generate an algebraic model, capable of being used with tdmore()
#'
#' @param predictFunction a function with at least arguments times and regimen. Any extra values are used as parameters.
#'
#' @return An algebraic model
#' @export
algebraic <- function(predictFunction) {
  argNames <- names( formals(predictFunction) )
  assertthat::assert_that(all(argNames[1:2] == c("times", "regimen")))
  pNames <- c()
  if(length(argNames) > 2) { pNames <- argNames[seq(3, length(argNames))] }

  structure(
    list(predictFunction = function(times, regimen, parameters) {
      do.call(predictFunction, c(list(times=times, regimen=regimen), parameters))
    },
    parameters=pNames
    ),
    class="algebraic"
  )
}

#' Predict for algebraic models
#'
#' We only support dosing into the default compartment, and only bolus doses
#'
#' @param model The model itself.
#' @param newdata dataframe with at least a column 'TIME' and other values. A prediction will be generated for each filled-in value.
#' @param regimen dataframe with column 'TIME' and adhering to standard NONMEM specifications otherwise (columns AMT, RATE, CMT)
#' @param parameters either a dataframe with column 'TIME' and a column for each covariate and parameter, or a named numeric vector
#'
#' @return
#' A data.frame similar to the newdata data frame, but with predicted values.
#'
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

  if(class(parameters) == "data.frame") {
    stop("Time-changing covariates not supported")
  } else {
    assertthat::assert_that(is.numeric(parameters))
    pNames <- names(parameters)
    params = parameters
  }
  assertthat::assert_that(all(pNames %in% model$parameters))
  assertthat::assert_that(all(model$parameters %in% pNames))

  # All arguments look good, let's prepare the simulation
  times <- newdata$TIME
  CONC <- model$predictFunction(times, regimen, params)

  # Only get the values we want
  result <- data.frame(TIME=times, CONC=CONC)
  result[, colnames(newdata)]
}


#' Build a tdmore object based on an algebraic model
#'
#' @param model The algebraic model
#' @param parameters Lits of parameters, should be NULL
#' @param add additive residual error, as stdev
#' @param prop proportional residual error, as stdev
#' @param exp exponential residual error, as stdev. The exponential error cannot be used in conjunction with the additive or proportional error
#'
#' @return a tdmore object, capable of estimating bayesian individual parameters
#' @export
tdmore.algebraic <- function(model, parameters=NULL, omega=NULL, add=0, prop=0, exp=0) {
  # Check that parameters + covariates together supplies the parameters needed for the model
  if(!is.null(parameters)) stop("Algebraic models can only work with their own parameters")

  assertthat::is.number(add)
  assertthat::is.number(prop)
  assertthat::is.number(exp)
  ## Exponential and add/prop are mutually exclusive
  if(exp != 0) assertthat::assert_that(add==0 & prop==0)
  if(add != 0 || prop != 0) assertthat::assert_that(exp == 0)

  structure(list(
    model=model,
    omega=omega,
    res_var=list(add=add, prop=prop, exp=exp),
    parameters=model$parameters
  ), class="tdmore")
}
