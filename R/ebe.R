#' Calculate the population log likelihood.
#'
#' @param estimate the current estimate of the parameters
#' @param omega the omega matrix of the model
#'
#' @return the population log likelihood
pop_ll <- function(estimate, omega) {
  sum( mvtnorm::dmvnorm(estimate, sigma=omega, log=TRUE) )
}

#' Calculate the prediction log likelihood.
#'
#' @param estimate the current estimate of the parameters
#' @param tdmore the tdmore object
#' @param observed data frame with at least a TIME column, and all observed data. The observed data will be compared to the model predictions.
#' If not specified, we estimate the population prediction
#' @param regimen data frame describing the treatment regimen.
#' @param covariates the model covariates
#'
#' @return the prediction log likelihood
pred_ll <- function(estimate, tdmore, observed, regimen, covariates) {
  if(is.null(observed) || nrow(observed) == 0) return(0)
  pred <- predict.tdmore(object=tdmore, newdata=observed, regimen=regimen, parameters=estimate, covariates=covariates)
  res <- residuals.tdmore(tdmore, observed, pred, log=TRUE)
  sum(res)
}

#' Calculate the log likelihood as the sum of the population likelihood and the prediction likelihood.
#'
#' @param estimate the current estimate of the parameters
#' @param tdmore the tdmore object
#' @param observed data frame with at least a TIME column, and all observed data. The observed data will be compared to the model predictions.
#' If not specified, we estimate the population prediction
#' @param regimen data frame describing the treatment regimen.
#' @param covariates the model covariates
#'
#' @return the log likelihood
ll <- function(estimate, tdmore, observed, regimen, covariates) {
  if(is.null(names(estimate))) names(estimate) <- tdmore$parameters #you should support named parameters as well!
  res <- pop_ll(estimate, tdmore$omega) + pred_ll(estimate, tdmore, observed, regimen, covariates)
  res
}


#' Calculate the Emperical Bayesian Estimated parameters that predict
#' a given dataset using the model
#'
#' @param tdmore a tdmore object
#' @param observed data frame with at least a TIME column, and all observed data. The observed data will be compared to the model predictions.
#' If not specified, we estimate the population prediction
#' @param regimen data frame describing the treatment regimen.
#' @param covariates the model covariates
#' @param p optional starting parameter for the MLE minimization
#' @param ... Extra parameters to pass to nlm
#'
#' @return A tdmorefit object
#' @importFrom stats nlm
#' @export
estimate <- function(tdmore, observed=NULL, regimen, covariates=NULL, p=NULL, ...) {
  assert_that(class(tdmore) == "tdmore")
  if(is.null(p)) p<-rep(0, length(tdmore$parameters))
  # First try to estimate at starting values, as a precaution
  ll(estimate=p,
     tdmore=tdmore,
     observed=observed,
     regimen=regimen,
     covariates=covariates)
  # Then do the full nlm
  pointEstimate <- stats::nlm(f=function(...) {
    tryCatch({
      -2*ll(...)
    }, error=function(e){ 999999 })
    },
                       p=p,
                       tdmore,
                       observed,
                       regimen,
                       covariates,
                       ...)
  res <- pointEstimate$estimate
  names(res) <- tdmore$parameters

  #OFIM <- pointEstimate$hessian ##!!! Wrong hessian!?
  # Observed fisher information matrix is -1*hessian(ll)
  OFIM <- numDeriv::hessian(ll, res, tdmore=tdmore, observed=observed, regimen=regimen, covariates=covariates)
  varcov <- solve(OFIM) #inverse of OFIM is an estimator of the asymptotic covariance matrix
  if(all(diag(varcov) <= 0)) varcov <- -1 * varcov
  dimnames(varcov) = list(tdmore$parameters, tdmore$parameters)
  ofv <- pointEstimate$minimum
  tdmorefit(tdmore, observed, regimen, covariates, ofv, res, varcov, nlmresult=pointEstimate)
}

#' Create a tdmorefit object manually
#'
#' @param tdmore the tdmore object
#' @param observed the observed data.frame, or NULL
#' @param regimen the treatment regimen data.frame
#' @param covariates the model covariates
#' @param ofv (optional) the OFV value
#' @param res the found parameter values, as a named vector, or NULL to use 0
#' @param varcov the found varcov matrix, or NULL to use a diagonal matrix
#' @param nlmresult the result of the non-linear minimization
#'
#' @return A tdmorefit object, manually created
#' @export
tdmorefit <- function(tdmore, observed=NULL, regimen, covariates=NULL, ofv=NA, res=NULL, varcov=NULL, nlmresult=NULL) {
  N <- length(tdmore$parameters)
  if(is.null(res)) {
    res <- rep(0, N) #population prediction
    names(res) <- tdmore$parameters
  }
  if(is.null(varcov)) {
    varcov <- diag(1, nrow=N, ncol=N)
    dimnames(varcov) = list(tdmore$parameters, tdmore$parameters)
  }

  structure(
    list(
      tdmore=tdmore,
      observed=observed,
      regimen=regimen,
      covariates=covariates,
      ofv=ofv,
      logLik=ofv / -2,
      res=res,
      varcov=varcov,
      nlmresult=nlmresult
    ), class=c("tdmorefit"))
}

#' Summary of an EBE fit
#'
#' @param object A tdmorefit object
#' @param ... Additional parameters ignored
#'
#' @return Nothing
#' @export
summary.tdmorefit <- function(object, ...) {
  tdmorefit <- object
  cat("Fit of model to ", nrow(tdmorefit$observed), " points of data\n")
  cat("LogLikelihood = ", logLik.tdmorefit(tdmorefit), "\n")
  cat("Point estimate: ", coef.tdmorefit(tdmorefit), "\n")
  cat("Variance/covariance matrix: ", vcov.tdmorefit(tdmorefit), "\n")
}

#' The obtained parameter values for the fit
#'
#' @param object A tdmorefit object
#' @param ... ignored
#'
#' @return Named numeric vector
#' @export
coef.tdmorefit <- function(object, ...) {object$res}

#' The variance-covariance matrix
#'
#' @param object A tdmorefit object
#' @param ... ignored
#'
#' @return Named numeric matrix
#' @export
vcov.tdmorefit <- function(object, ...) {object$varcov}

#' Calculates the confidence interval on the obtained parameter values
#'
#' @param object A tdmorefit object
#' @param ... ignored
#'
#' @return A matrix (or vector) with columns giving
#' lower and upper confidence limits for each parameter.
#' These will be labelled as (1-level)/2
#' and 1 - (1-level)/2 in % (by default 2.5% and 97.5%).
#' @export
confint.tdmorefit <- function(object, ...) {
  stats::confint.default(object)
}

#' Get the predicted (fitted) data
#'
#' @param object A tdmorefit object
#' @param ... ignored
#'
#' @return A data.frame
#' @export
fitted.tdmorefit <- function(object, ...) {
  predict.tdmore(object=object$tdmore, newdata=object$observed, regimen=object$regimen, parameters=object$res, covariates=object$covariate)
}

#' Predict new data using the current model
#'
#' @param object A tdmorefit object
#' @param newdata
#' A data.frame with new data and the columns to predict,
#' or a numeric vector to specify times, and predict all model output
#' or NULL to interpolate between 0 and the maximum known times
#' @param regimen Treatment regimen
#' @param parameters Set parameters. If missing, or if only part of the parameters are specified,
#' the other parameters are taken from the tdmorefit object
#' @param covariates the model covariates, named vector, or data.frame with column 'TIME', and at least TIME 0
#' @param se.fit TRUE to provide a confidence interval on the prediction, adding columns xxx.median, xxx.upper and xxx.lower
#' @param level The confidence interval, or NA to return all mc.maxpts results
#' @param mc.maxpts Maximum number of points to sample in Monte Carlo simulation
#' @param ip.maxpts Maximum number of points to interpolate if newdata is not given
#' @param .progress see plyr::ddply
#' @param .parallel see plyr::ddply
#' @param ... ignored
#'
#' @return A data.frame
#' @export
#'
#' @importFrom stats coef vcov median quantile
predict.tdmorefit <- function(object, newdata=NULL, regimen=NULL, parameters=NULL, covariates=NULL, se.fit=FALSE, level=0.95, mc.maxpts=100, ip.maxpts=100, .progress="none", .parallel=FALSE, ...) {
  tdmorefit <- object
  if(is.null(regimen)) regimen=tdmorefit$regimen
  if(is.null(newdata)) {
    tmax = regimen$TIME + regimen$II * regimen$ADDL
    tmax = max(0, tmax, tdmorefit$observed$TIME)
    newdata = seq(0, tmax, length.out=ip.maxpts)
  }
  if(is.numeric(newdata)) {
    # keep as numeric
  } else {
    assert_that("TIME" %in% colnames(newdata))
  }
  pars <- coef(tdmorefit)
  if(!is.null(parameters)) {
    pars[names(parameters)] <- parameters ## set pars from argument
  }
  if(is.null(covariates)) covariates <- tdmorefit$covariates
  ipred <- predict.tdmore(object=tdmorefit$tdmore, newdata=newdata, regimen=regimen, parameters=pars, covariates=covariates)
  if(se.fit) {
    oNames <- names(newdata)
    oNames <- oNames[oNames != "TIME"]
    pNames <- names(coef(tdmorefit))
    mc <- as.data.frame( mnormt::rmnorm(mc.maxpts, mean=coef(tdmorefit), varcov=vcov(tdmorefit)) )
    mc$sample <- 1:mc.maxpts
    for(i in names(parameters)) mc[, i] <- parameters[i]
    fittedMC <- plyr::ddply(mc, 1, function(row) {
      res <- row[pNames]
      pred <- predict.tdmore(object=tdmorefit$tdmore, newdata=newdata, regimen=regimen, parameters=unlist(res), covariates=covariates)
      cbind(row, pred)
      #pred$sample <- row$sample
      #pred
    }, .progress=.progress, .parallel=.parallel)
    if(is.na(level)) { #user requested full dataframe without summary
      return(fittedMC)
    }
    a <- (1-level)/2
    plyr::ddply(fittedMC, "TIME", function(x) {
      result <- list(TIME = x$TIME[1])
      for(i in oNames) {
        result[i] <- ipred[ipred$TIME==x$TIME[1], i]
        result[paste0(i, ".median")] <- median(x[,i])
        result[paste0(i, ".lower")] <- quantile(x[,i], probs=a)
        result[paste0(i, ".upper")] <- quantile(x[,i], probs=1-a)
      }
      unlist(result)
    })
  } else {
    ipred
  }
}

#' Get the log-likelihood of the predicted values,
#' or try to calculate it
#'
#' @param object A tdmorefit object
#' @param ... ignored
#'
#' @return A numeric value
#' @export
logLik.tdmorefit <- function(object, ...) {
  res <- object$logLik
  if(is.finite(res))
    return(res)
  # Damn, it was not previously calculated
  estimate <- coef(object)
  tdmore <- formula(object)
  observed <- model.frame(object)
  regimen <- object$regimen
  covariates <- object$covariates
  ll(estimate=estimate,
     tdmore=tdmore,
     observed=observed,
     regimen=regimen,
     covariates=covariates)
}

#' Get an overview of the log-likelihood for varying parameter values
#'
#' @param fitted A tdmorefit object
#' @param maxpts Maximum number of points per parameter
#' @param limits Limits to explore (numeric vector of form c(min, max))
#' @param .progress See plyr::ddply
#' @param ... ignored
#'
#' @return A data.frame with each parameter value tested, and an additional `logLik` column with the log-likelihood for each parameter combination
#' @export
profile.tdmorefit <- function(fitted, maxpts = 100, limits=c(-1.96, 1.96), .progress="none", ...) {
  tdmorefit <- fitted
  model <- tdmorefit$tdmore
  x <- seq(limits[1], limits[2], length.out=maxpts)
  x <- expand.grid(rep(list(x), length(model$parameters)))
  colnames(x) <- model$parameters
  result <- plyr::adply(x, 1, function(estimate) {
    c(logLik=ll(as.numeric(estimate), model, tdmorefit$observed, tdmorefit$regimen, tdmorefit$covariates))
  }, .progress=.progress)
  result
}

#' Get the model used to provide this fit
#'
#' @param x A tdmorefit object
#' @param ... ignored
#'
#' @return The original tdmore object
#' @export
formula.tdmorefit <- function(x, ...) {x$tdmore}

#' Get the observed values used to provide this model fit.
#'
#' @param formula A tdmorefit object
#' @param se if TRUE, add a column xx.upper and xx.lower with the lower and upper bounds of confidence, based on the residual error model
#' @param level Confidence interval
#' @param ... ignored
#'
#' @return A data.frame, similar to the original
#' @export
model.frame.tdmorefit <- function(formula, se=FALSE, level=0.95, ...) {
  tdmorefit <- formula
  if(is.null(tdmorefit$observed)) return(NULL)
  result <- model.frame.tdmore(tdmorefit$tdmore, tdmorefit$observed, se=se, level=level)
  result
}
