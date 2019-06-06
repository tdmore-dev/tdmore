#'
#' Add an iteration column 'ITER' to the observed dataframe.
#' Observations will be grouped by occasions.
#'
#' @param regimen regimen
#' @param observed observed
#'
#' @return modified observed dataframe, with an 'ITER' column
#'
addIterationColumn <- function(regimen, observed) {
  observedMpc <- observed
  observedMpc$OCC <- 1
  for (obsIndex in seq_len(nrow(observed))) {
    obsTime <- observed$TIME[obsIndex]
    observedMpc$OCC[obsIndex] <- regimen %>% filter(TIME < obsTime) %>% pull(OCC) %>% last()
  }
  iterationRow <- !duplicated(observedMpc$OCC)
  observedMpc[iterationRow, "ITER"] <- seq_len(sum(iterationRow))
  observedMpc$ITER <- zoo::na.locf(observedMpc$ITER)
  return(observedMpc %>% dplyr::select(-one_of("OCC")))
}

#'
#' Select most appropriate covariates.
#'
#' @param t time
#' @param covariates time-varying covariates
#' @param thetaNames theta names to be removed
#' @return a vector with the covariates
#'
selectBestCovariates <- function(t, covariates, thetaNames) {
  selectedRow <- utils::tail(which((covariates$TIME <= t)==TRUE), n=1)
  retValue <- unlist(covariates[selectedRow, ] %>% dplyr::select(-one_of(thetaNames, "TIME")))
  return(retValue)
}

#'
#' MPC estimate method.
#'
#' @inheritParams estimate
#' @return a tdmorefit
#' @export
#'
estimate.tdmore_mpc <- function(object, observed=NULL, regimen=NULL, covariates=NULL, par=NULL, fix=NULL, method="L-BFGS-B", se.fit=TRUE, lower=NULL, upper=NULL, multistart=F, control=list(), data=NULL, ...) {
  if (is.numeric(covariates)) {
    covariates <- as.data.frame(as.list(c(TIME=0, covariates)))
  }

  # Init
  theta <- object$mpc_theta
  thetaNames <- names(theta)
  thetaDf <- as.data.frame(as.list(theta)) # Thetas (MPC parameters) specified to tdmore via the covariates dataframe
  thetaDf$TIME <- 0
  covariatesMpc <- cbind(thetaDf, t(selectBestCovariates(t=0, covariates, thetaNames)))
  fix <- c()

  # Special case: no observed data
  if (is.null(observed) || nrow(observed) == 0) {
    ipred <- estimate.tdmore(object, regimen=regimen, covariates=covariatesMpc, se.fit=TRUE, control=control, ...)
    maxIter <- 0
  # Normal case: observed data present
  } else {
    observedMpc <- addIterationColumn(regimen, observed)
    maxIter <- max(observedMpc$ITER)
  }

  for (iterIndex in seq_len(maxIter + 1)) {
    if (iterIndex > 1) {
      firstDoseJustAfter <- regimen %>% filter(TIME >= currentObsTime) %>% dplyr::pull(TIME) %>% dplyr::first()
      if (is.na(firstDoseJustAfter)) {
        break
      }
      previousEbe <- predict(ipred, newdata=c(firstDoseJustAfter))[, paste0(thetaNames, object$mpc_suffix)]
      names(previousEbe) <- paste0(thetaNames)
      previousEbe$TIME <- firstDoseJustAfter
      covariates_tmp <- cbind(previousEbe, t(selectBestCovariates(t=firstDoseJustAfter, covariates, thetaNames)))
      covariatesMpc <- dplyr::bind_rows(covariatesMpc, covariates_tmp)
      fix <- coef(ipred)
    }
    if (iterIndex <= maxIter) {
      currentObsTime <- observedMpc %>% dplyr::filter(ITER==iterIndex) %>% dplyr::pull(TIME) %>% dplyr::last()
      ipred <- estimate.tdmore(object, regimen=regimen %>% dplyr::filter(TIME < currentObsTime),
                        observed=observedMpc %>% dplyr::filter(ITER==iterIndex) %>% dplyr::select(-one_of("ITER")),
                        covariates=covariatesMpc, fix=fix, se.fit=se.fit, control=control, ...)
    } else {
      # Last extra iteration copies final covariates dataframe (with all MPC parameters) to ipred
      ipred$covariates <- covariatesMpc
    }
  }
  return(ipred)
}

#' MPC is a generic function to make a tdmore model compatible with MPC.
#'
#' @param x a model object (type not defined)
#' @param theta named vector with the initial values
#' @param suffix a suffix to add to the theta names
#' @param ... extra arguments
#' @return an object of class tdmore_mpc
#' @export
mpc <- function(x, theta, suffix, ...) {
  UseMethod("mpc")
}

#'
#' Build a 'tdmore_mpc' object based on a 'tdmore' object.
#'
#' @inheritParams mpc
#' @return an object of class tdmore_mpc
#' @export
mpc.tdmore <- function(x, theta, suffix, ...) {
  x$mpc_theta <- theta
  x$mpc_suffix <- suffix
  class(x) <- append("tdmore_mpc", "tdmore")
  return(x)
}
