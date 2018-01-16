## Library of functions
pk1cptivbolusVCL <- function(THETA=list(V=10, CL=5), OMEGA=list(V=0.20, CL=0.30)) {
  function(times, regimen, ETA_V, ETA_CL) {
    assert_that(all(colnames(regimen) == c("TIME", "AMT")))
    V = THETA$V * exp(OMEGA$V * ETA_V)
    CL = THETA$CL * exp(OMEGA$CL * ETA_CL)
    k = CL / V
    t = times

    CONC <- rep(0, length(times))
    for(i in seq(1, nrow(regimen))) {
      tD = regimen$TIME[i]
      D = regimen$AMT[i]
      CONC <- CONC + ifelse(t >= tD, D / V * exp(-k*(t-tD)), 0)
    }
    CONC
  }
}

algebraic <- function(predictFunction) {
  argNames <- names( formals(predictFunction) )
  assert_that(all(argNames[1:2] == c("times", "regimen")))
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
#' @param observed dataframe with at least a column 'TIME' and other values. A prediction will be generated for each filled-in value.
#' @param regimen dataframe with column 'TIME' and adhering to standard NONMEM specifications otherwise (columns AMT, RATE, CMT)
#' @param parameters either a dataframe with column 'TIME' and a column for each covariate and parameter, or a named numeric vector
#'
#' @return
#' A data.frame similar to the observed data frame, but with predicted values.
#'
model_predict.algebraic <- function(model, observed, regimen=data.frame(TIME=c()), parameters=c()) {
  # Verify arguments are good
  assert_that("data.frame" %in% class(observed))
  assert_that("TIME" %in% colnames(observed))
  oNames <- colnames(observed)
  oNames <- oNames[oNames != "TIME"]
  assert_that(length(oNames) > 0)
  assert_that(all(oNames %in% c("CONC")))

  assert_that("data.frame" %in% class(regimen))
  assert_that(all(c("TIME", "AMT") %in% colnames(regimen)))
  assert_that(all(colnames(regimen) %in% c("TIME", "AMT", "RATE")))

  if(class(parameters) == "data.frame") {
    stop("Time-changing covariates not supported")
  } else {
    assert_that(is.numeric(parameters))
    pNames <- names(parameters)
    params = parameters
  }
  assert_that(all(pNames %in% model$parameters))
  assert_that(all(model$parameters %in% pNames))

  # All arguments look good, let's prepare the simulation
  times <- observed$TIME
  CONC <- model$predictFunction(times, regimen, params)

  # Only get the values we want
  result <- data.frame(TIME=times, CONC=CONC)
  result[, colnames(observed)]
}


tdmore.algebraic <- function(model, parameters=NULL, add=0, prop=0, exp=0) {
  # Check that parameters + covariates together supplies the parameters needed for the model
  if(!is.null(parameters)) stop("Algebraic models can only work with their own parameters")

  assertthat::is.number(add)
  assertthat::is.number(prop)
  assertthat::is.number(exp)
  ## Exponential and add/prop are mutually exclusive
  if(exp != 0) assert_that(add==0 & prop==0)
  if(add != 0 || prop != 0) assert_that(exp == 0)

  structure(list(
    model=model,
    res_var=list(add=add, prop=prop, exp=exp),
    parameters=model$parameters
  ), class="tdmore")
}
