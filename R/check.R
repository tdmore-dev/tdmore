checkTdmore <- function(tdmore) {
  assert_that(is.tdmore(tdmore))
  tdmore <- checkModel(tdmore)
  tdmore <- checkOmegaMatrix(tdmore)
  tdmore <- checkNamesConsistency(tdmore)
  tdmore <- checkErrorModel(tdmore)
  tdmore <- checkIOV(tdmore)
  return(tdmore)
}

#' @importFrom utils str
checkOmegaMatrix <- function(tdmore) {
  omega <- tdmore$omega
  parameters <- tdmore$parameters

  if(is.null(omega)) omega = diag(rep(1, length(parameters)))
  if(!is.numeric(omega)) stop("`omega' should be a numeric vector or matrix")
  if(is.null(dim(omega))) omega <- vectorToDiagonalMatrix(omega)

  assertthat::are_equal(nrow(omega), length(parameters))
  assertthat::are_equal(ncol(omega), length(parameters))
  assertthat::assert_that(all( eigen(omega)$values >= 0), msg = "Omega matrix is not positive semi-definite")
  assert_that(all(diag(omega)!=0), msg = "Omega's can't have a zero value")

  colnames <- colnames(omega)
  rownames <- rownames(omega)

  if (is.null(colnames) || is.null(rownames)) {
    # If omega provided without any information regarding the column/row names
    colnames(omega) <- parameters
    rownames(omega) <- parameters
  } else {
    assert_that(all(parameters %in% colnames),
                msg = paste0("Inconsistent omega column (", utils::str(colnames) , ") and parameter names (", utils::str(parameters) ,")"))
    assert_that(all(parameters %in% rownames),
                msg = paste0("Inconsistent omega row (", utils::str(rownames), ") and parameter names (", utils::str(parameters), ")"))

    # Reorder omega column and row names
    omega <- omega[, parameters, drop=FALSE]
    omega <- omega[parameters, , drop=FALSE]
  }

  ## Double check that omega is symmetric
  assert_that( isSymmetric(omega), msg="Omega matrix not symmetric")

  tdmore$omega <- omega
  return(tdmore)
}

checkErrorModel <- function(tdmore) {
  for (errorModel in tdmore$res_var) {
    add <- errorModel$add
    prop <- errorModel$prop

    assertthat::is.number(add)
    assertthat::is.number(prop)

    if (inherits(tdmore$model, "RxODE")) {
      modVars <- tdmore$model$get.modelVars()
      assert_that(errorModel$var %in% c(modVars$lhs, modVars$state),
                  msg = paste0("Error model variable `", errorModel$var, "' is not predicted by the model"))
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

  if (inherits(model, "RxODE")) {
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
  if(is.data.frame(covariates)) {
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

checkIOV <- function(tdmore) {
  iov <- unlist(tdmore$iov)
  if (!is.null(iov)) {
    if (is.character(iov)) {
      parameters <- tdmore$parameters
      condition <- iov %in% parameters
      assert_that(all(condition), msg=paste("IOV term(s)", paste(iov[!condition], collapse=",")  ,"not defined in model"))
      tdmore$iov <- parameters[(parameters %in% iov)] # Reorder IOV params in iov attribute

      ## Check that no correlation exists in OMEGA
      iovI <- parameters %in% iov
      assert_that( all(tdmore$omega[ iovI, !iovI ] == 0), msg="No correlation between IOV and non-IOV term can exist in OMEGA." )
    } else {
      stop("IOV should be a character or character array")
    }
  }
  return(tdmore)
}

checkRegimen <- function(regimen, iov) {
  assert_that(is.data.frame(regimen))
  assert_that(all(c("TIME", "AMT") %in% colnames(regimen)))
  assert_that(all(colnames(regimen) %in% c("TIME", "AMT", "RATE", "DURATION", "CMT", "II", "ADDL", "SS", "OCC")))

  # XXX: Not a strict requirement! II + SS is also possible...
  #if ("II" %in% colnames(regimen) || "ADDL" %in% colnames(regimen))
  #  assert_that(all(c("II", "ADDL") %in% colnames(regimen)))

  if ("OCC" %in% colnames(regimen)) {
    if (is.null(iov)) {
      stop("No IOV exists in model")
    } else {
      uniqueOcc <- unique(regimen$OCC)
      assert_that(all(uniqueOcc == 1:max(regimen$OCC)), msg="Occasions mispecified in regimen")
    }
  }
  if (!("OCC" %in% colnames(regimen)) && !is.null(iov)) {
    stop("No occasion column (OCC) in regimen")
  }

}

checkTimes <- function(times) {
  assert_that(is.numeric(times))
  assert_that(!is.unsorted(times))
}
