checkTdmore <- function(tdmore) {
  assert_that("tdmore" %in% class(tdmore))
  tdmore <- checkModel(tdmore)
  tdmore <- checkOmegaMatrix(tdmore)
  tdmore <- checkNamesConsistency(tdmore)
  tdmore <- checkErrorModel(tdmore)
  return(tdmore)
}

checkOmegaMatrix <- function(tdmore) {
  omega <- tdmore$omega
  parameters <- tdmore$parameters

  if(is.null(omega)) omega = diag(rep(1, length(parameters)))

  assertthat::are_equal(nrow(omega), length(parameters))
  assertthat::are_equal(ncol(omega), length(parameters))
  assertthat::assert_that(all( eigen(omega)$values >= 0), msg = "Omega matrix is not positive semi-definite")
  assert_that(sum(diag(tdmore$omega)==0) == 0, msg = "Omega's can't have a zero value")

  colnames <- colnames(omega)
  rownames <- rownames(omega)

  if (is.null(colnames) || is.null(rownames)) {
    # If omega provided without any information regarding the column/row names
    colnames(omega) <- parameters
    rownames(omega) <- parameters
  } else {
    assert_that(all(parameters %in% colnames), msg = "Inconsistent omega column and parameter names")
    assert_that(all(parameters %in% rownames), msg = "Inconsistent omega row and parameter names")

    # Reorder omega column and row names
    omega <- omega[, parameters, drop=FALSE]
    omega <- omega[parameters, , drop=FALSE]
  }

  tdmore$omega <- omega
  return(tdmore)
}

checkErrorModel <- function(tdmore) {
  for (errorModel in tdmore$res_var) {
    add <- errorModel$add
    prop <- errorModel$prop
    exp <- errorModel$exp

    assertthat::is.number(add)
    assertthat::is.number(prop)
    assertthat::is.number(exp)

    if (exp != 0) assert_that(add==0 & prop==0, msg = "Exponential and add/prop are not mutually exclusive")
    if (add != 0 || prop != 0) assert_that(exp == 0, msg = "Exponential and add/prop are not mutually exclusive")

    if (c("RxODE") %in% class(tdmore$model)) {
      modVars <- tdmore$model$get.modelVars()
      assert_that(errorModel$var %in% c(modVars$lhs, modVars$state), msg = "Unknown variable defined in error model")
    }
  }
  return(tdmore)
}

checkNamesConsistency <- function(tdmore) {
  parameters <- tdmore$parameters
  omegaNames <- names(diag(tdmore$omega))
  assert_that(all(parameters == omegaNames), msg = "Inconsistent omega and parameter names")
  return(tdmore)
}

checkModel <- function(tdmore) {
  parameters <- tdmore$parameters
  model <- tdmore$model
  covariates <- tdmore$covariates

  if (c("RxODE") %in% class(model)) {
    # Check that parameters + covariates together supplies the parameters needed for the model
    modVars <- model$get.modelVars()
    if(is.null(parameters)) parameters <- modVars$params
    if(is.null(covariates)) covariates <- modVars$params[ ! modVars$params %in% parameters ] #the rest
    assert_that( all.equal( sort(c(parameters, covariates)), sort(modVars$params)), msg = "Inconsistent model, some covariates seem to be missing")
  }

  tdmore$parameters <- parameters
  tdmore$covariates <- covariates
  return(tdmore)
}

checkCovariates <- function(tdmore, covariates) {
  if("data.frame" %in% class(covariates)) {
    cNames <- colnames(covariates)
    cNames <- cNames[cNames != "TIME"]
  } else if (is.numeric(covariates)){
    cNames <- names(covariates)
  } else if (length(covariates) == 0) {
    cNames <- c()
  } else {
    stop("Covariates in wrong format")
  }
  condition <- tdmore$covariates %in% cNames
  assert_that(all(condition), msg=paste("Value for covariate(s)", paste(tdmore$covariates[!condition], collapse=",")  ,"missing"))
}
