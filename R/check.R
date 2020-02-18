checkTdmore <- function(tdmore) {
  stopifnot(is.tdmore(tdmore))
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

  stopifnot(nrow(omega) == length(parameters))
  stopifnot(ncol(omega) == length(parameters))
  if(!all( eigen(omega)$values >= 0)) stop("Omega matrix is not positive semi-definite")
  if(!all(diag(omega)!=0)) stop("Omega's can't have a zero value")

  colnames <- colnames(omega)
  rownames <- rownames(omega)

  if (is.null(colnames) || is.null(rownames)) {
    # If omega provided without any information regarding the column/row names
    colnames(omega) <- parameters
    rownames(omega) <- parameters
  } else {
    if(!all(parameters %in% colnames))
                stop(paste0("Inconsistent omega column (", utils::str(colnames) , ") and parameter names (", utils::str(parameters) ,")"))
    if(!all(parameters %in% rownames))
                stop(paste0("Inconsistent omega row (", utils::str(rownames), ") and parameter names (", utils::str(parameters), ")"))

    # Reorder omega column and row names
    omega <- omega[, parameters, drop=FALSE]
    omega <- omega[parameters, , drop=FALSE]
  }

  ## Double check that omega is symmetric
  if( !isSymmetric(omega)) stop("Omega matrix not symmetric")

  tdmore$omega <- omega
  return(tdmore)
}

checkErrorModel <- function(tdmore) {
  for (errorModel in tdmore$res_var) {
    add <- errorModel$add
    prop <- errorModel$prop

    stopifnot(is.numeric(add))
    stopifnot(is.numeric(prop))

    if (inherits(tdmore$model, "RxODE")) {
      modVars <- tdmore$model$get.modelVars()
      if(! errorModel$var %in% c(modVars$lhs, modVars$state))
                  stop(paste0("Error model variable `", errorModel$var, "' is not predicted by the model"))
    }
  }
  return(tdmore)
}

checkNamesConsistency <- function(tdmore) {
  parameters <- tdmore$parameters
  omegaNames <- names(diag(tdmore$omega))
  if(!all(parameters == omegaNames)) stop("Inconsistent omega and parameter names")
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
    if(!isTRUE(all.equal( sort(c(parameters, covariates)), sort(modVars$params)))) stop("Inconsistent model, some covariates seem to be missing")
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
  if(!all(condition)) stop(paste("Value for covariate(s)", paste(tdmore$covariates[!condition], collapse=",")  ,"missing"))
}

checkIOV <- function(tdmore) {
  iov <- unlist(tdmore$iov)
  if (!is.null(iov)) {
    if (is.character(iov)) {
      parameters <- tdmore$parameters
      condition <- iov %in% parameters
      if(!all(condition)) stop(paste("IOV term(s)", paste(iov[!condition], collapse=",")  ,"not defined in model"))
      tdmore$iov <- parameters[(parameters %in% iov)] # Reorder IOV params in iov attribute

      ## Check that no correlation exists in OMEGA
      iovI <- parameters %in% iov
      if( !all(tdmore$omega[ iovI, !iovI ] == 0)) stop("No correlation between IOV and non-IOV term can exist in OMEGA." )
    } else {
      stop("IOV should be a character or character array")
    }
  }
  return(tdmore)
}

checkRegimen <- function(regimen, iov) {
  stopifnot(is.data.frame(regimen))
  if(!all(c("TIME", "AMT") %in% colnames(regimen)))
              stop("A required column is missing in the treatment regimen")

  # XXX: Not a strict requirement! II + SS is also possible...
  #if ("II" %in% colnames(regimen) || "ADDL" %in% colnames(regimen))
  #  stopifnot(all(c("II", "ADDL") %in% colnames(regimen)))

  if ("OCC" %in% colnames(regimen)) {
    if (is.null(iov)) {
      shouldWarn <- getOption("tdmore.warnIov", NA)
      if(is.na(shouldWarn)) {
        warning("No IOV exists in model\nThis warning will appear once during the runtime of this application. To enable consistently, set options(tdmore.warnIov)")
        options(tdmore.warnIov=FALSE)
      } else if (shouldWarn) {
        warning("No IOV exists in model")
      } else {
        # do not warn
      }
      #XXX: this is not an issue...
    } else {
      uniqueOcc <- unique(regimen$OCC)
      if(!all(uniqueOcc == 1:max(regimen$OCC)))
        stop("Occasions mispecified in regimen")
    }
  }
  if (!("OCC" %in% colnames(regimen)) && !is.null(iov) && length(iov) > 0) {
    stop("No occasion column (OCC) in regimen")
  }
}

checkTimes <- function(times) {
  stopifnot(is.numeric(times))
  stopifnot(!is.unsorted(times))
}
