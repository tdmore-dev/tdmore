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
  observedMpc <- addIterationColumn(regimen, observed)
  theta <- object$theta
  thetaNames <- names(theta)
  thetaDf <- as.data.frame(as.list(theta)) # Thetas (MPC parameters) specified to tdmore via the covariates dataframe
  thetaDf$TIME <- 0
  maxIter <- max(observedMpc$ITER)

  for (iterIndex in seq_len(maxIter + 1)) {
    if (iterIndex == 1) {
      covariatesMpc <- cbind(thetaDf, t(selectBestCovariates(t=0, covariates, thetaNames)))
      fix <- c()
    } else {
      prevObsTime <- currentObsTime
      previousEbe <- predict(ipred, newdata=c(prevObsTime))[, paste0(thetaNames, object$suffix)]
      names(previousEbe) <- paste0(thetaNames)
      previousEbe$TIME <- prevObsTime
      covariates_tmp <- cbind(previousEbe, t(selectBestCovariates(t=prevObsTime, covariates, thetaNames)))
      covariatesMpc <- dplyr::bind_rows(covariatesMpc, covariates_tmp)
      fix <- coef(ipred)
    }
    if (iterIndex <= maxIter) {
      currentObsTime <- observedMpc %>% dplyr::filter(ITER==iterIndex) %>% dplyr::pull(TIME) %>% dplyr::last()
      ipred <- estimate(object$model, regimen=regimen %>% dplyr::filter(TIME < currentObsTime),
                        observed=observedMpc %>% dplyr::filter(ITER==iterIndex) %>% dplyr::select(-one_of("ITER")),
                        covariates=covariatesMpc, fix=fix, se.fit=FALSE, ...)
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
  structure(list(
    model=x,
    theta=theta,
    suffix=suffix
  ), class="tdmore_mpc")
}
