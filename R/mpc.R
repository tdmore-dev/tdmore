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
  i <- findInterval(observed$TIME, regimen$TIME, left.open=T)
  OCC <- regimen$OCC[i]
  observed$ITER <- as.numeric( factor(OCC) ) #count from 1 to N
  observed
}

#'
#' MPC predict method.
#'
#' @inheritParams predict.tdmore
#' @return a data.frame with all observed values at the given time points
#' @export
#'
predict.tdmore_mpc <- function(object, newdata, regimen=NULL, parameters=NULL, covariates=NULL, se=FALSE, level=0.95, ...) {
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(as.list(c(TIME=0, covariates)))
  }

  # Retrieve initial thetas
  theta <- object$mpc_theta
  thetaNames <- names(theta)

  # Make distinction between covariates
  allCovariates <- object$covariates
  trueCovariates <- allCovariates[!(allCovariates %in% thetaNames)]
  assert_that(all(trueCovariates %in% colnames(covariates)), msg="Missing covariates")

  # Add missing initial thetas if not there
  missingThetas <- thetaNames[!(thetaNames %in% colnames(covariates))]
  for (missingTheta in missingThetas) {
    covariates[, missingTheta] <- theta[[missingTheta]]
  }

  return(predict.tdmore(object=object, newdata=newdata, regimen=regimen, parameters=parameters, covariates=covariates, se=se, level=level, ...))
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

  # Retrieve initial thetas
  theta <- object$mpc_theta
  thetaNames <- names(theta)

  # Make distinction between covariates
  allCovariates <- object$covariates
  trueCovariates <- allCovariates[!(allCovariates %in% thetaNames)]
  assert_that(sum(thetaNames %in% colnames(covariates))==0, msg="You shouldn't give thetas in covariates")
  assert_that(all(trueCovariates %in% colnames(covariates)), msg="Some covariates are missing")

  # Initialisation
  ebeCovariates <- as.data.frame(as.list(theta)) # Thetas (MPC parameters) specified to tdmore via the covariates dataframe
  ebeCovariates$TIME <- 0
  fix <- c()

  # Special case: no observed data
  if (is.null(observed) || nrow(observed) == 0) {
    # Note that se.fit=T => varcov can then be used to plot BSV in population
    ipred <- estimate.default(object, regimen=regimen, covariates=mergeCovariates(ebeCovariates, covariates), se.fit=TRUE, control=control, ...)
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
      ebeCovariates <- dplyr::bind_rows(ebeCovariates, previousEbe)
      fix <- coef(ipred)
    }
    if (iterIndex <= maxIter) {
      currentObsTime <- observedMpc %>% dplyr::filter(.data$ITER==iterIndex) %>% dplyr::pull(.data$TIME) %>% dplyr::last()
      ipred <- estimate.default(object, regimen=regimen %>% dplyr::filter(.data$TIME < currentObsTime),
                        observed=observedMpc %>% dplyr::filter(.data$ITER==iterIndex) %>% dplyr::select(-dplyr::one_of("ITER")),
                        covariates=mergeCovariates(ebeCovariates, covariates), fix=fix, se.fit=se.fit, control=control, ...)
    } else {
      # Last extra iteration copies final covariates dataframe (with all MPC parameters) to ipred
      ipred$covariates <- mergeCovariates(ebeCovariates, covariates)
    }
  }

  # Put the original data back in the `observed` and `regimen` slot
  ipred$observed <- observed
  ipred$regimen <- regimen

  return(ipred)
}

#' Merge covariates.
#'
#' @param ebeCovariates EBE covariates
#' @param covariates 'true' time-varying covariates
#' @return an updated covariates dataframe
#'
mergeCovariates <- function(ebeCovariates, covariates) {
  join <- dplyr::full_join(ebeCovariates, covariates, by="TIME")
  join <- join[order(join$TIME), ]
  return(zoo::na.locf(join))
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
  if( !setequal(x$parameters, x$iov) ) stop("For MPC, all parameters should be IOV")
  return(x)
}
