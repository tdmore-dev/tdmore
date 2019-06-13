#'
#' Add an iteration column 'ITER' to the observed dataframe.
#' The iteration column determines in which MPC iteration each observation
#' will be included.
#'
#' This is calculated by matching the corresponding OCC from the regimen (time-wise).
#'
#' @param regimen regimen
#' @param observed observed
#'
#' @return modified observed dataframe, with an 'ITER' column
#'
addIterationColumn <- function(regimen, observed) {
  i <- findInterval(observed$TIME, regimen$TIME)
  OCC <- regimen$OCC[i]
  observed$ITER <- as.numeric( factor(OCC) ) #count from 1 to N
  observed
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
  retValue <- unlist(covariates[selectedRow,,drop=F] %>% dplyr::select(-dplyr::one_of(thetaNames, "TIME")))
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
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(as.list(c(TIME=0, covariates)))
  }

  # Init
  theta <- object$mpc_theta
  thetaNames <- names(theta)
  thetaDf <- as.data.frame(as.list(theta)) # Thetas (MPC parameters) specified to tdmore via the covariates dataframe
  thetaDf$TIME <- 0
  selectedCovs <- selectBestCovariates(t=0, covariates, thetaNames)
  covariatesMpc <- if(is.null(selectedCovs)) {thetaDf} else {cbind(thetaDf, t(selectedCovs))}
  fix <- c()

  # Special case: no observed data
  if (is.null(observed) || nrow(observed) == 0) {
    # Note that se.fit=T => varcov can then be used to plot BSV in population
    ipred <- estimate.tdmore(object, regimen=regimen, covariates=covariatesMpc, se.fit=TRUE, control=control, ...)
    maxIter <- 0
  # Normal case: observed data present
  } else {
    observedMpc <- addIterationColumn(regimen, observed)
    maxIter <- max(observedMpc$ITER)
  }

  for (iterIndex in seq_len(maxIter + 1)) {
    if (iterIndex > 1) {
      firstDoseJustAfter <- regimen$TIME[ match(TRUE, regimen$TIME >= currentObsTime) ]
      if (is.na(firstDoseJustAfter)) {
        break
      }
      previousEbe <- predict(ipred, newdata=c(firstDoseJustAfter))[, paste0(thetaNames, object$mpc_suffix)]
      names(previousEbe) <- paste0(thetaNames)
      previousEbe$TIME <- firstDoseJustAfter
      selectedCovs <- selectBestCovariates(t=firstDoseJustAfter, covariates, thetaNames)
      covariates_tmp <- if(is.null(selectedCovs)) {previousEbe} else {cbind(previousEbe, t(selectedCovs))}
      covariatesMpc <- dplyr::bind_rows(covariatesMpc, covariates_tmp)
      fix <- coef(ipred)
    }
    if (iterIndex <= maxIter) {
      currentObsTime <- observedMpc %>% dplyr::filter(.data$ITER==iterIndex) %>% dplyr::pull(.data$TIME) %>% dplyr::last()
      ipred <- estimate.tdmore(object, regimen=regimen %>% dplyr::filter(.data$TIME < currentObsTime),
                        observed=observedMpc %>% dplyr::filter(.data$ITER==iterIndex) %>% dplyr::select(-dplyr::one_of("ITER")),
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