#' Calculate the population log likelihood.
#'
#' @param par the current estimate of the parameters
#' @param omega the omega matrix
#' @param tdmore the tdmore object
#' @param observed the observed data, not used in the population log likelihood
#' @param regimen data frame describing the treatment regimen
#' @param covariates the model covariates, not used in the population log likelihood
#'
#' @return the population log likelihood
pop_ll <- function(par, omega, tdmore, observed, regimen, covariates) {
  sum <- sum( mvtnorm::dmvnorm(par, sigma=omega, log=TRUE) )
  return(sum)
}

#' Calculate the prediction log likelihood.
#'
#' @param par the current estimate of the parameters
#' @param omega the omega matrix
#' @param tdmore the tdmore object
#' @param observed data frame with at least a TIME column, and all observed data.
#' The data.frame can be empty, or could contain only a TIME column.
#' The observed data will be compared to the model predictions.
#' If not specified, we estimate the population prediction
#' @param regimen data frame describing the treatment regimen
#' @param covariates the model covariates
#'
#' @importFrom stats residuals
#'
#' @return the prediction log likelihood
pred_ll <- function(par, omega, tdmore, observed, regimen, covariates) {
  pred <- predict(object=tdmore, newdata=observed, regimen=regimen, parameters=par, covariates=covariates)
  res <- residuals(tdmore, observed, pred, log=TRUE)
  return(sum(res))
}

#' Calculate the log likelihood as the sum of the population likelihood and the prediction likelihood.
#'
#' @param par the current estimate of the parameters
#' @param omega the omega matrix
#' @param tdmore the tdmore object
#' @param observed data frame with at least a TIME column, and all observed data. The observed data will be compared to the model predictions.
#' If not specified, we estimate the population prediction
#' @param regimen data frame describing the treatment regimen.
#' @param covariates the model covariates
#'
#' @return the log likelihood
ll <- function(par, omega, tdmore, observed, regimen, covariates) {
  names(par) <- getParameterNames(tdmore, regimen)
  pop_ll <- pop_ll(par, omega, tdmore, observed, regimen, covariates)
  pred_ll <- pred_ll(par, omega, tdmore, observed, regimen, covariates)
  return(pop_ll + pred_ll)
}

#' Calculate the Empirical Bayesian Estimate parameters that predict
#' a given dataset using the model
#'
#' @param object a tdmore, tdmore_set or tdmore_mixture object
#' @param observed data frame with at least a TIME column, and all observed data. The observed data will be compared to the model predictions.
#' If not specified, we estimate the population prediction
#' @param regimen data frame describing the treatment regimen.
#' @param covariates the model covariates
#' @param par optional starting parameter for the MLE minimization
#' @param method the optimisation method, by default, method "L-BFGS-B" is used
#' @param se.fit calculate the variance-covariance matrix for the fit
#' Putting this to FALSE can reduce computation time.
#' @param lower the lower bounds of the parameters, if null, -5 * sqrt(diag(model$omega)) is used
#' @param upper the upper bounds of the parameters, if null, +5 * sqrt(diag(model$omega)) is used
#' @param ... extra parameters to pass to the optim function
#' A good example is to specify `control=list(trace=1, REPORT=10, factr=1e13)` to improve performance.
#' This will `trace` the estimation progress
#' and `report` every 10 iterations. It limits the precision of the `L-BFGS-B` method to 3 significant digits (`factr` is multiplied with the machine precision, `1e12 * 1e-16 = 1e-4`).
#'
#' Incidentally, 3 significant digits \href{http://holford.fmhs.auckland.ac.nz/research/sigdig}{is the default value for NONMEM}.
#' Without the above option, `estimate` will calculate with a relative precision of `1e-8`.
#'
#' @return A tdmorefit object
#' @importFrom stats optim
#' @export
estimate <- function(object, observed=NULL, regimen, covariates=NULL, par=NULL, method="L-BFGS-B", se.fit=TRUE, lower=NULL, upper=NULL, ...) {

  if (inherits(object, "tdmore_set")) {
    # Either a tdmore or tdmore_mixture
    object <- findFirstCompatibleModel(object, covariates)
  }
  if (inherits(object, "tdmore")) {
    tdmore <- object
  }
  else if (inherits(object, "tdmore_mixture")) {
    return(estimateMixtureModel(object, observed=observed, regimen=regimen, covariates=covariates, par=par, method=method, lower=lower, upper=upper, ...))
  }
  else {
    stop("Object not an instance of 'tdmore', 'tdmore_set' or 'tdmore_mixture")
  }
  observed <- model.frame(tdmore, data=observed) #ensure "observed" in right format for estimation

  # Processing IOV
  iov <- tdmore$iov
  occasions <- getMaxOccasion(regimen)
  parNames <- getParameterNames(tdmore, regimen)

  # Setting optim initial conditions
  par <- processParameters(par, tdmore, regimen)
  omega <- expandOmega(tdmore, occasions)
  lower <- processParameters(lower, tdmore, regimen, -5 * sqrt(diag(omega)))
  upper <- processParameters(upper, tdmore, regimen, +5 * sqrt(diag(omega)))

  # First try to estimate at starting values, as a precaution
  value <- ll(par=par, omega=omega, tdmore=tdmore, observed=observed, regimen=regimen, covariates=covariates)
  if(!is.finite(value))
    stop("Log-likelihood is ", value, " at starting values `par`. Cannot start optimization routine.")

  # Function to optimise
  fn <- function(par, ...) {
    tryCatch({
      -2 * ll(par = par, omega = omega, ...)
    }, error = function(e) {
      999999
    })
  }

  # Then do the full optimisation
  pointEstimate <- stats::optim(
    par = par,
    fn = fn,
    method = method,
    lower = lower,
    upper = upper,
    hessian = se.fit, #only calculate the hessian if we want the vcov
    tdmore = tdmore,
    observed = observed,
    regimen = regimen,
    covariates = covariates,
    ...
  )
  res <- pointEstimate$par
  names(res) <- parNames

  # Observed fisher information matrix = -hessian(ll)
  if(se.fit) {
    OFIM <- pointEstimate$hessian * 1/2
    varcov <- solve(OFIM) #inverse of OFIM is an estimator of the asymptotic covariance matrix
  } else {
    varcov <- diag(0, nrow=length(parNames)) #vcov not calculated, so assumed '0'
  }
  dimnames(varcov) = list(parNames, parNames)
  ofv <- pointEstimate$value
  tdmorefit(tdmore, observed, regimen, covariates, ofv, res, varcov, nlmresult=pointEstimate, call=match.call())
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
#' @param call a informative string that gives the arguments that were used in estimate()
#'
#' @return A tdmorefit object, manually created
#' @export
tdmorefit <- function(tdmore, observed=NULL, regimen, covariates=NULL, ofv=NA, res=NULL, varcov=NULL, nlmresult=NULL, call=NULL) {
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
      logLik=ofv/-2,
      res=res,
      varcov=varcov,
      nlmresult=nlmresult,
      call=call
    ), class=c("tdmorefit"))
}

#' Print an EBE fit.
#'
#' @param x a tdmorefit object
#' @param ... ignored
#'
#' @export
print.tdmorefit <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("Coef:\n")
  print(x$res)
}


#' Summary of a tdmorefit object.
#'
#' @param object a tdmorefit object
#' @param ... additional parameters ignored
#' @importFrom stats logLik
#'
#' @export
summary.tdmorefit <- function(object, ...) {
  x <- object
  cat("Call:\n")
  print(x$call)
  cat("Coef:\n")
  print(x$res)
  cat("\n")
  cat(paste0(
    "OFV: ",
    signif(logLik(object, "ll") * -2, digits=3),
    " (pop=",
    signif(logLik(object, "pop") * -2, digits=3),
    ", pred=",
    signif(logLik(object, "pred") * -2, digits=3),
    ")\n"
  ))
  cat("\n")
  cat("Observations:\n")
  print(x$observed, row.names=F)
  cat("\n")
  cat("Regimen:\n")
  print(x$regimen, row.names=F)
  cat("\n")
  cat("Covariates:", if(is.null(x$covariates)) "/\n" else "\n")
  if(!is.null(x$covariates)) {
    print(x$covariates, row.names=F)
  }
  cat("\n")
  cat("Coefficients:\n")
  coefDf <- data.frame(name=names(x$res), value=signif(x$res, digits=3), se=signif(sqrt(diag(x$varcov)), digits=3))
  coefDf$value.low <- signif(coefDf$value - 1.96*coefDf$se, digits=3)
  coefDf$value.up <- signif(coefDf$value + 1.96*coefDf$se, digits=3)
  coefDf$ci <- paste0('(', coefDf$value.low, ', ', coefDf$value.up, ')')
  coefDf$value.low <- NULL
  coefDf$value.up <- NULL
  colnames(coefDf) <- c('name', 'value', 'se', '(95% CI)')
  coefDf <- coefDf[,c(1,2,4,3)]
  print(coefDf, row.names=F)
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
#' @param .progress see plyr::ddply
#' @param .parallel see plyr::ddply
#' @param ... ignored
#'
#' @return A data.frame
#' @export
#'
#' @importFrom stats coef vcov median quantile
predict.tdmorefit <- function(object, newdata=NULL, regimen=NULL, parameters=NULL, covariates=NULL, se.fit=FALSE, level=0.95, mc.maxpts=100, .progress="none", .parallel=FALSE, ...) {
  tdmorefit <- object
  if(is.null(regimen)) regimen=tdmorefit$regimen
  if(is.null(newdata)) newdata <- model.frame(tdmorefit)

  pars <- processParameters(parameters, tdmorefit$tdmore, regimen, defaultValues=coef(tdmorefit))
  if(is.null(covariates)) covariates <- tdmorefit$covariates

  ipred <- predict(object=tdmorefit$tdmore, newdata=newdata, regimen=regimen, parameters=pars, covariates=covariates)
  if(se.fit) {
    mc <- as.data.frame( mnormt::rmnorm(mc.maxpts, mean=coef(tdmorefit), varcov=vcov(tdmorefit)) )
    colnames(mc) <- names(pars)
    mc <- cbind( sample=1:mc.maxpts, mc ) #make sure 'sample' is first column
    uniqueColnames <- make.unique(colnames(mc)) # needed for dplyr to have unique colnames

    for(i in names(parameters)) mc[, i] <- parameters[i] # TODO: adapt because of same names
    fittedMC <- plyr::ddply(mc, 1, function(row) {
      res <- unlist(row[-1]) # Remove 'sample'
      names(res) <- names(pars)
      pred <- predict(object=tdmorefit$tdmore, newdata=newdata, regimen=regimen, parameters=res, covariates=covariates)
      colnames(row) <- uniqueColnames
      resArray <- cbind(row, pred)
      resArray
    }, .progress=.progress, .parallel=.parallel)

    colnames(fittedMC)[seq_len(length(uniqueColnames))] <- colnames(mc)

    if(is.na(level)) { #user requested full dataframe without summary
      return(fittedMC)
    }
    oNames <- getPredictOutputNames(newdata, colnames(fittedMC), names(coef(tdmorefit)))
    return(summariseFittedMC(fittedMC, ipred, level, oNames))

  } else {
    ipred
  }
}

#' Return the output names of the fitted MC matrix.
#'
#' @param newdata the newdata data frame
#' @param columnNames current column names
#' @param pNames parameters names
#'
#' @return a vector with all the output names
getPredictOutputNames <- function(newdata, columnNames, pNames) {
  oNames <- names(newdata)
  oNames <- oNames[oNames != "TIME"]

  # By default, if newdata is numeric vector (oNames is null), all outputs from RxODE are returned
  if(is.null(oNames)) {
    oNames <- columnNames[!(columnNames %in% c("time", "TIME", "fit", "sample", pNames))]
  }
  return(oNames)
}

#' Summarise fitted MC matrix over all the samples.
#' TIME column must be present.
#'
#' @param fittedMC fitted MC matrix (result of predict.tdmorefit)
#' @param ipred result of predict.tdmorefit
#' @param level ignored
#' @param oNames outputs to summarise
#'
#' @return a summarised data frame
summariseFittedMC <- function(fittedMC, ipred, level, oNames) {
  a <- (1-level)/2
  retValue <- plyr::ddply(fittedMC, "TIME", function(x) {
    result <- list(TIME = x$TIME[1])
    for(i in oNames) {
      result[i] <- ipred[ipred$TIME==x$TIME[1], i]
      result[paste0(i, ".median")] <- median(x[,i])
      result[paste0(i, ".lower")] <- quantile(x[,i], probs=a)
      result[paste0(i, ".upper")] <- quantile(x[,i], probs=1-a)
    }
    unlist(result)
  })
  return(retValue)
}

#' Get the log-likelihood of the predicted values.
#'
#' @param object A tdmorefit object
#' @param type log-lokelihood function type, 3 possible values: 'pop', 'pred' or 'll' (= pop + pred)
#' @param ... ignored
#'
#' @return A numeric value
#' @export
#' @importFrom stats formula model.frame
logLik.tdmorefit <- function(object, type=c('ll', 'pop', 'pred'), ...) {
  estimate <- coef(object)
  tdmore <- formula(object)
  observed <- model.frame(object)
  regimen <- object$regimen
  covariates <- object$covariates
  fun <- getLikelihoodFun(type)
  fun(par=estimate,
      omega=expandOmega(tdmore, getMaxOccasion(regimen)),
      tdmore=tdmore,
      observed=observed,
      regimen=regimen,
      covariates=covariates)
}

#' Get an overview of the log-likelihood for varying parameter values.
#' If arguments `fix` and `limits` are not set for a specific parameter, individual eta's will be generated to cover 95 percent of the population distribution.
#'
#' @param fitted A tdmorefit object
#' @param fix Which parameters to fix? Named vector of the parameters that are fixed and should not be profiled
#' @param maxpts Maximum number of points per parameter
#' @param limits limits to explore (numeric vector of form c(min, max) or specific limits per parameter in the form list(ETA1=c(min,max), ETA2=c(min,max), etc))
#' @param type log-lokelihood function type, 3 possible values: 'pop', 'pred' or 'll' (= pop + pred)
#' @param .progress See plyr::ddply
#' @param ... ignored
#'
#' @return a tdmore profile object. It namely contains a data.frame with each parameter value tested, and an additional `logLik` column with the log-likelihood for each parameter combination.
#' @export
profile.tdmorefit <- function(fitted, fix=NULL, maxpts = 50, limits=NULL, type=c('ll', 'pop', 'pred'), .progress="none", ...) {
  tdmorefit <- fitted
  model <- tdmorefit$tdmore
  omegas <- diag(tdmorefit$tdmore$omega)
  profiledParameters <- model$parameters
  if(!is.null(fix)) profiledParameters <- profiledParameters[ !(profiledParameters %in% names(fix)) ]

  limitsAsList <- !is.null(limits) && is.list(limits)
  if(limitsAsList) assert_that(all(names(limits) %in% profiledParameters))

  limitsAsNumeric <- !is.null(limits) && is.numeric(limits)
  if(limitsAsNumeric) assert_that(length(limits)==2)

  fun <- getLikelihoodFun(type)

  list <- lapply(
    profiledParameters,
    FUN = function(profiledParameter) {
        if (limitsAsNumeric) {
          return(seq(limits[1], limits[2], length.out=maxpts))

        } else if (limitsAsList){
          parameterRange <- unlist(limits[profiledParameter])
          if(!is.null(parameterRange)) {
            cat()
            return(seq(parameterRange[1], parameterRange[2], length.out=maxpts))
          }
        }
        # Default case, take omega value to generate a 95% range of possible eta's
        omega <- omegas[profiledParameter]
        return(seq(-1.96 * sqrt(omega), 1.96 * sqrt(omega), length.out = maxpts))
    }
  )
  grid <- expand.grid(list)
  colnames(grid) <- profiledParameters
  for(i in names(fix)) grid[,i] <- fix[i]
  grid <- grid[, model$parameters, drop=FALSE] # Reorder the columns

  profile <- plyr::adply(grid, 1, function(estimate) {
    eta <- as.numeric(estimate)
    names(eta) <- model$parameters
    omega <- expandOmega(model, getMaxOccasion(tdmorefit$regimen))
    c(logLik=fun(eta, omega, model, tdmorefit$observed, tdmorefit$regimen, tdmorefit$covariates))
  }, .progress=.progress)

  return(structure(
    list(
      profile = profile,
      profiledParameters = profiledParameters,
      tdmorefit = tdmorefit
    ),
    class = c("tdmoreprofile")
  ))
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
#' @param data Data.frame to append and modify. If NULL, uses the observed values.
#' @param se if TRUE, add a column xx.upper and xx.lower with the lower and upper bounds of confidence, based on the residual error model
#' @param level Confidence interval to use for `se`
#' @param ... ignored
#'
#' @return
#' A data.frame, similar to the one used to estimate this `tdmorefit` object.
#' If `se` was specified, then a column xx.upper and xx.lower with the
#' lower and upper confidence interval (based on the residual error model) is added.
#'
#' @export
model.frame.tdmorefit <- function(formula, data=NULL, se=FALSE, level=0.95, ...) {
  if(is.null(data)) data <- formula$observed
  result <- model.frame.tdmore(formula=formula$tdmore, data=data, se=se, level=level)
  result
}

#' Get the right likelihood function for the specified type.
#'
#' @param type likelihood type, ll, pop or pred
#'
#' @return the correct tdmore function
getLikelihoodFun <- function(type) {
  type <- match.arg(type, c('ll', 'pop', 'pred'))
  if(type == "ll") fun <- ll
  else if (type == "pop") fun <- pop_ll
  else if (type == "pred") fun <- pred_ll
  else stop("Unknown log-likelihood function")
}

#' Check if the given argument is of class `tdmorefit`
#' @param x An object to test
#' @keywords internal
#' @export
is.tdmorefit <- function(a) {inherits(a, "tdmorefit")}

#' Expand OMEGA matrix by adding variance-covariance information regarding all IOV terms.
#' Note that the column and row names of the returned matrix are strictly identical to the ones returned by getParameterNames().
#'
#' @param tdmore the tdmore object
#' @param occasions how many occasions
#' @return the omega matrix with the duplicated IOV terms at the end
#' @export
expandOmega <- function(tdmore, occasions) {
  iov <- tdmore$iov
  omega <- tdmore$omega

  if(is.null(iov)) {
    retValue <- omega
  } else {
    indexes <- seq_len(ncol(omega))
    iovIndexes <- indexes[(colnames(omega) %in% iov)]
    noIovIndexes <- indexes[!(colnames(omega) %in% iov)]
    mat_tmp <- omega[c(noIovIndexes, iovIndexes), c(noIovIndexes, iovIndexes)]

    for(occasion in seq_len(occasions - 1)) {
      indexesToCopy <- (ncol(mat_tmp)-(length(iov) - 1)):ncol(mat_tmp)
      # Copy vertical IOV columns at the end
      matV <- mat_tmp[,indexesToCopy]
      mat_tmp <- cbind(mat_tmp, matV)

      # Copy horizontal IOV columns at the end
      matH <- mat_tmp[indexesToCopy,]
      mat_tmp <- rbind(mat_tmp, matH)

      # Remove initial IOV correlations
      for(sourceIndex in indexesToCopy) {
        destIndex <- sourceIndex + length(iov)
        mat_tmp[sourceIndex, destIndex] <- 0
        mat_tmp[destIndex, sourceIndex] <- 0
      }
    }
    retValue <- mat_tmp
  }
  return(retValue)
}
