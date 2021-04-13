#' Calculate the population log likelihood. Fixed values are ignored.
#'
#' @inheritParams ll
#' @return the population log likelihood
#' @engine
pop_ll <- function(par, omega, fix, tdmore, observed, regimen, covariates, isChol=FALSE) {
  stopifnot(all.equal(rownames(omega), names(par)) )

  if(isTRUE( getOption("tdmore.mvtnorm", default=FALSE)) && requireNamespace("mvnfast") ) {
    pdf <- mvnfast::dmvn( X=par, mu=par*0, sigma=omega, log=TRUE, isChol=isChol )
  } else {
    if(isChol) {
      omega <- solve(chol2inv(omega))
    }
    pdf <- mvtnorm::dmvnorm(x=par, mean=par*0, sigma=omega, log=TRUE)
  }

  sum <- sum( pdf )
  return(sum)
}

#' Calculate the prediction log likelihood. Fixed values are appended to the provided parameter estimate
#'
#' @inheritParams ll
#' @return the prediction log likelihood
#' @engine
pred_ll <- function(par, omega, fix, tdmore, observed, regimen, covariates) {
  stopifnot(all.equal(rownames(omega), names(par)) )
  pred <- stats::predict(object=tdmore, newdata=observed, regimen=regimen, parameters=c(fix$values,par), covariates=covariates)
  ll <- 0
  for(res_var in tdmore$res_var) {
    i <- res_var$var
    if(! i %in% colnames(pred)) next #not in observed values, ignore
    ipred <- pred[, i, drop=TRUE]
    obs <- observed[, i, drop=TRUE]
    thisLL <- res_var$ll(ipred, obs)
    ll <- ll + sum(thisLL, na.rm=TRUE )
  }
  if(ll==0) {
    vars <- sapply(tdmore$res_var, function(x){x$var})
    if(nrow(observed)>0 && !any(vars %in% colnames(observed))) {
      warning("Observed has rows, but none contributed to pred_ll. Are you definitely sure an output column is present? ", paste(vars))
    }
  }
  return(ll)
}

#' Calculate the log likelihood as the sum of the population likelihood and the prediction likelihood.
#'
#' @param par a named numeric vector with parameter values; fixed parameter values are not included
#' @param omega the omega matrix; fixed parameter values are not included
#' @param fix a named list with value 'indexes' specifying the fixed indexes, and value 'values' specifying a named numeric vector with the actual fixed values.
#' @param tdmore the tdmore object
#' @param observed data frame with at least a TIME column, and all observed data. The observed data will be compared to the model predictions.
#' If not specified, we estimate the population prediction
#' @param regimen data frame describing the treatment regimen.
#' @param covariates the model covariates
#' @param isChol was omega specified as the Cholesky decomposition?
#'
#' @return the log likelihood
#' @engine
ll <- function(par, omega, fix, tdmore, observed, regimen, covariates, isChol=FALSE) {
  stopifnot(all.equal(rownames(omega), names(par)) )

  parNames <- getParameterNames(tdmore, regimen)
  if(!is.null(fix$indexes)) {
    parNames <- parNames[-fix$indexes] #names without the fixIndexes
    names(par) <- parNames
  }
  pop_ll <- pop_ll(par, omega, fix, tdmore, observed, regimen, covariates, isChol=isChol)
  if(is.null(observed) || nrow(observed) == 0) {
    pred_ll <- 0
  } else {
    pred_ll <- pred_ll(par, omega, fix, tdmore, observed, regimen, covariates)
  }
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
#' @param fix named numeric vector with the fixed parameters. Names can repeat to fix parameters in multiple occasions.
#' @param method the optimisation method, by default, method "L-BFGS-B" is used
#' Can also be specified as a list, in which case all methods specified will be tried and the one finding the estimate with highest log-likelihood will be returned.
#' @param se.fit calculate the variance-covariance matrix for the fit
#' Putting this to FALSE can reduce computation time.
#' @param lower the lower bounds of the parameters, if null, -5 * sqrt(diag(model$omega)) is used
#' @param upper the upper bounds of the parameters, if null, +5 * sqrt(diag(model$omega)) is used
#' @param multistart if TRUE, perform optimization using many different starting conditions.
#' A combination of -sd, +sd for each non-IOV parameter is explored.
#' This avoids returning a local optimum, and strengthens the belief
#' that the found optimum is a global optimum.
#'
#' Alternatively, a named vector with the requested perturbations to parameters
#' can be provided: \code{estimate(multistart=c(eta_Tlag=0.5))}. This will evaluate
#' the provided starting value for eta_Tlag +- 0.5, without also perturbing the other initial values.
#'
#' WARNING: This may lead to extremely long computation times. It is recommended to use a control argument to
#' reduce accuracy and report estimation progress.
#'
#' @param control control options for `estimate`
#' If trace > 0, the estimation function also reports errors when calculating log-likelihood.
#' @param .mpc if FALSE, assume normal EBE estimation
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
#' @engine
estimate <- function(object, observed, regimen, covariates, par, fix,
                             method="L-BFGS-B", se.fit=TRUE,
                             lower=NA, upper=NA,
                             multistart=F,
                             control=list(trace=getOption("tdmore.trace", default=(!testthat::is_testing()&interactive())*1), REPORT=10, factr=1e13),
                             ...,
                             .progress=NULL, .mpc=TRUE) {
  if(missing(observed)) observed <- NULL
  if(missing(regimen)) regimen <- NULL
  if(missing(covariates)) covariates <- NULL
  if(missing(par)) par <- NULL
  if(missing(fix)) fix <- NULL
  if(is.null(control$trace)) control$trace <- 0

  if(.mpc && inherits(object, "tdmore_mpc")) {
    return( estimate_tdmore_mpc(
      object=object, observed=observed, regimen=regimen, covariates=covariates, par=par, fix=fix,
      method=method, se.fit=se.fit,
      lower=lower, upper=upper,
      multistart=multistart,
      control=control,
      ...,
      .progress=.progress
    ))
  }

  if(length(method) > 1) {
    ## Multiple methods specified
    ## Try them one at a time
    ipred <- NULL
    bestLL <- -Inf
    call <- match.call()
    for(m in method) {
      call$method <- m
      #thisIpred <- eval(call)
      thisIpred <- estimate(
        object, observed, regimen, covariates, par, fix, method=m, se.fit, lower, upper, multistart, control, ...
      )
      thisLL <- logLik(thisIpred)
      if(thisLL > bestLL) {
        if(control$trace > 0) cat(m," found better fit (LL=", thisLL, ")\n")
        bestLL <- thisLL
        ipred <- thisIpred
      } else {
        if(control$trace > 0) cat(m," not better (LL=", thisLL, ")\n")
      }
    }
    return(ipred)
  }

  if (is.tdmore(object)) {
    tdmore <- object
  } else if (inherits(object, "tdmore_mixture")) {
    return(estimateMixtureModel(object, observed=observed, regimen=regimen, covariates=covariates, par=par, fix=fix, method=method, lower=lower, upper=upper, ...))
  } else if (is.tdmorefit(object)) {
    ## Re-estimate with different parameters/options

    # if using an MPC object, remove the theta names from the stored covariates first
    covariates <- covariates %||% object$covariates

    return(
      estimate(object$tdmore,
               observed %||% object$observed,
               regimen %||% object$regimen,
               covariates,
               par=par, fix=fix, method=method, se.fit=se.fit, lower=lower, upper=upper, multistart=multistart, control=control, ...)
      )
  } else {
    #will never happen
    stop("Object not an instance of 'tdmore', 'tdmore_set' or 'tdmore_mixture") #nocov
  }

  observed <- model.frame(tdmore, data=observed) #ensure "observed" in right format for estimation

  # Processing IOV
  iov <- tdmore$iov
  occasions <- getMaxOccasion(regimen)
  parNames <- getParameterNames(tdmore, regimen)

  # Setting optim initial conditions
  par <- processParameters(par, tdmore, regimen)
  omega <- expandOmega(tdmore, occasions)
  if(!is.null(lower)) lower <- processParameters(lower, tdmore, regimen, -5 * sqrt(diag(omega)))
  if(!is.null(upper)) upper <- processParameters(upper, tdmore, regimen, +5 * sqrt(diag(omega)))

  # Prepare the tdmore cache
  cTdmore <- tdmore
  if(is.null(cTdmore$cache)) {
    cTdmore$cache <- model_prepare(tdmore$model, times=observed$TIME, regimen=regimen, parameters=par, covariates=covariates, iov=tdmore$iov, extraArguments=tdmore$extraArguments)
  } else {
    ## TODO: Check if cache is still valid?
  }

  # Process fix vector
  fixIndexes <- getFixedParametersIndexes(parNames=names(par), fix=fix)
  allIndexes <- seq_len(length(par))
  parIndexes <- allIndexes[! allIndexes %in% fixIndexes]
  par <- par[parIndexes]
  omega <- omega[parIndexes, parIndexes, drop=FALSE]
  if(!is.null(lower)) lower <- lower[parIndexes]
  if(!is.null(upper)) upper <- upper[parIndexes]
  fix <- list(indexes=fixIndexes, values=fix)
  updatedParNames <- parNames[parIndexes]

  # First try to estimate at starting values, as a precaution
  value <- ll(par=par, omega=omega, fix=fix, tdmore=cTdmore, observed=observed, regimen=regimen, covariates=covariates)
  if(!is.finite(value))
    stop("Log-likelihood is ", value, " at starting values `par`. Cannot start optimization routine.")

  # Function to optimise
  fn <- function(par, omega, fix, tdmore, observed, regimen, covariates, isChol) {
    tryCatch({
      val <- -2 * ll(par=par, omega=omega, fix=fix, tdmore=tdmore, observed=observed, regimen=regimen, covariates=covariates, isChol=isChol)
      #cat("LL calculated for ", str(par), ": ", val, "\n")
      val
    }, error = function(e) {
      if(control$trace > 0) print(e)
      999999
    })
  }

  omegaChol <- Matrix::chol(omega)

  # Then do the full optimisation
  arg <- list(
    par = par,
    fn = fn,
    method = method,
    lower = lower,
    upper = upper,
    hessian = se.fit, #only calculate the hessian if we want the vcov
    tdmore = cTdmore,
    observed = observed,
    regimen = regimen,
    covariates = covariates,
    isChol=TRUE,
    omega=omegaChol,
    fix=fix,
    control=control,
    ...
  )
  if(nrow(observed)==0) {
    control <- c(list(maxit=0), control) #no iterations, just use default!
  }

  if(!isFALSE(multistart)) {
    if(isTRUE(multistart)) {
      multistart <- sqrt( diag(omega)[ names(arg$par) ] )
      names(multistart) <- arg$par
    } else {
      stopifnot(all(names(multistart) %in% names(arg$par)))
    }

    parmat <- do.call( expand.grid, lapply(names(multistart), function(x){arg$par[x] + c(-1, 1) * multistart[x]}) )
    colnames(parmat) <- names(multistart)
    # The other ones do not change

    for(i in setdiff(names(arg$par), colnames(parmat))) {
      parmat[,i] <- arg$par[i]
    }
    parmat <- parmat[, names(arg$par), drop=F] #ensure same order
    multiArg <- arg
    if(ncol(parmat)==1) {
      # multistart drops the parameter names if only a single column is provided
      # see https://github.com/cran/optimr/issues/1
      multiArg$fn <- function(par, ...) {
        names(par) <- colnames(parmat)
        fn(par, ...)
      }
    }
    multiArg$parmat <- parmat
    multiArg$par <- NULL
    multiArg$hessian <- FALSE
    pointEstimates <- do.call(optimr::multistart, multiArg)
    pointEstimate <- pointEstimates[ which.min(pointEstimates$value), ,drop=FALSE]
    arg$par <- unlist(pointEstimate[ 1, seq_along(par), drop=FALSE] ) ## use found 'global' optimum
  }
  pointEstimate <- do.call(optimr::optimr, arg)
  res <- pointEstimate$par
  names(res) <- updatedParNames

  if(anyNA(res)) stop("Method ", method, " resulted in NA coefficients... Not returning result.")

  # Observed fisher information matrix = hessian(ll)
  # See https://stats.stackexchange.com/questions/68080/basic-question-about-fisher-information-matrix-and-relationship-to-hessian-and-s
  # OFIM is therefore -H/2, as we were optimizing -2*LL
  # And therefore the variance-covariance matrix is OFIM^-1
  if(se.fit) {
    ## Recalculate hessian?
    H <- pointEstimate$hessian / 2 ## we were optimizing -2*LL. The OFIM is the hessian for LL
    OFIM <- H
    varcov <- solve(OFIM) #inverse of OFIM is an estimator of the asymptotic covariance matrix
    if(varcov[1,1] < 0) varcov <- -1 * varcov
    tryCatch( {chol(varcov)}, #try the cholesky decomposition, to double check varcov is truly semi-positive definite!
              error=function(e) {
                varcov <<- NULL #set varcov to NULL, signifying that MCMC methods should be used to sample uncertainty
              })
  } else {
    varcov <- diag(.Machine$double.eps, nrow=length(updatedParNames)) #very small value, to keep matrix semi-definite
  }

  # Re-add fixed parameters in res and omega
  updatedRes <- rep(0, length(parNames))
  names(updatedRes) <- parNames
  updatedRes[parIndexes] <- res
  for(name in unique(names(fix$values))) {
    indexesInFix <- which(names(fix$values)==name)
    indexesInUpdatedRes <- which(names(updatedRes) == name)
    updatedRes[indexesInUpdatedRes[1:length(indexesInFix)]] <- fix$values[indexesInFix]
  }
  if(is.null(fix$indexes)) {
    updatedVarcov <- varcov # nothing to do
  } else {
    if(is.null(varcov)) {
      updatedVarcov <- NULL
    } else {
      updatedVarcov <- matrix(0, nrow=length(parNames), ncol=length(parNames))
      updatedVarcov[-fixIndexes, -fixIndexes] <- varcov
    }
  }

  if(!is.null(updatedVarcov)) dimnames(updatedVarcov) = list(parNames, parNames)
  ofv <- pointEstimate$value
  tdmorefit(tdmore, observed, regimen, covariates, ofv, updatedRes, updatedVarcov, fix$values, nlmresult=pointEstimate, call=match.call())
}

#' Get the indexes of the fixed parameters.
#'
#' @param parNames parameter names in the right order, possibly with repeated values
#' @param fix named numeric vector with the fixed parameters (user input)
#'
#' @keywords internal
#'
#' @return the fixed parameters indexes
#' @noRd
getFixedParametersIndexes <- function(parNames, fix) {
  if(length(fix) == 0) {
    return(c())
  }
  ## Boolean vector to store which parameters are fixed
  isFixed <- logical(length=length(parNames))

  for(name in unique(parNames)) {
    ## Does this name appear in the 'fix' parameters, and how often?
    N <- sum( names(fix) == name )
    if(N == 0) next #not fixed
    i <- which(parNames == name) # Which 'isFixed' should be set to TRUE?
    if(N > length(i)) stop("Inconsistent fix argument")
    i <- i[1:N] # only as many as we have 'fix' values for
    isFixed[i] <- TRUE
  }
  if(all(isFixed)) stop("Parameters cannot be all fixed")

  which(isFixed) # return indexes
}

#' Create a tdmorefit object manually
#'
#' @param tdmore the tdmore object
#' @param observed the observed data.frame, or NULL
#' @param regimen the treatment regimen data.frame
#' @param covariates the model covariates
#' @param ofv (optional) the OFV value
#' @param res the found parameter values, as a named vector, or NULL to use 0.
#' @param varcov the found varcov matrix
#' @param fix the fixed parameters
#' @param nlmresult the result of the non-linear minimization
#' @param call a informative string that gives the arguments that were used in estimate()
#'
#' @return A tdmorefit object, manually created
#' @engine
tdmorefit <- function(tdmore, observed=NULL, regimen, covariates=NULL, ofv=NA, res=NULL, varcov=NULL, fix=NULL, nlmresult=NULL, call=NULL) {
  N <- length(tdmore$parameters)
  if(is.null(res)) {
    res <- rep(0, N) #population prediction
    names(res) <- tdmore$parameters
  }

  ## Check parameter input
  extraNames <- setdiff(names(res), tdmore$parameters)
  if( !is.null(extraNames) && length(extraNames) > 0 ) stop("Superfluous names in parameter specification: ", paste(extraNames, collapse=", "), " does not appear in model parameter specification")
  missingNames <- setdiff(tdmore$parameters, names(res))
  if( !is.null(missingNames) && length(missingNames) > 0 ) stop("Missing parameters: ", paste(missingNames, collapse=", "))

  # if(is.null(varcov)) {
  #   varcov <- diag(.Machine$double.eps, nrow=N, ncol=N)
  #   dimnames(varcov) = list(tdmore$parameters, tdmore$parameters)
  # }

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
      fix=fix,
      nlmresult=nlmresult,
      call=call
    ), class=c("tdmorefit"))
}

#' @export
print.tdmorefit <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("Coef:\n")
  print(x$res)
}


#' @importFrom stats logLik
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

#' @export
coef.tdmorefit <- function(object, ...) {object$res}

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
#' @noRd
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
#' @noRd
#' @export
fitted.tdmorefit <- function(object, ...) {
  stats::predict(object=object$tdmore, newdata=object$observed, regimen=object$regimen, parameters=object$res, covariates=object$covariate)
}

#' Generate parameter samples using uncertainty around the fit.
#'
#' If a variance-covariance matrix could be generated during estimation using asymptotic , a simple multivariate normal sampling is used.
#' If the variance-covariance matrix is NULL, we use Metropolis-Hastings sampling.
#'
#' @param tdmorefit a tdmorefit object
#' @param fix named numeric vector with the fixed parameters
#' @param mc.maxpts number of Monte-Carlo samples
#' @param .progress Allows to specify a plyr-like progress object
#' A plyr progress object is a list with 3 function definitions: `init(N)`, `step()` and `term()`.
#' This can also be specified as a boolean. TRUE uses the default dplyr progress_estimated.
#'
#' @return the Monte-Carlo matrix, first column is always the 'sample' column
#' @engine
sampleMC <- function(tdmorefit, regimen=NULL, fix=tdmorefit$fix, mc.maxpts=100, .progress=interactive()) {
  mcmcAvailable <- requireNamespace("MCMCpack", quietly = TRUE)
  mnormtAvailable <- requireNamespace("mnormt", quietly = TRUE)
  if(!mnormtAvailable) {
    warning("The mnormt package is required to sample uncertainty from the FIM. As it is not installed, we will try to use Metropolis-Hastings instead.")
  }
  if(is.null(tdmorefit$varcov) || !mnormtAvailable) {
    if (!mcmcAvailable) stop("The MCMCpack must be installed to sample uncertainty using Metropolis-Hastings")
    return(
      sampleMC_metrop(tdmorefit, regimen=regimen, fix, mc.maxpts, .progress=.progress)
    )
  }
  sampleMC_norm(tdmorefit, regimen=regimen, fix, mc.maxpts)
}

#'
#' Generate parameter samples using a muilti-variate normal distribution around the variance-covariance matrix
#'
#' @param tdmorefit a tdmorefit object
#' @param fix named numeric vector with the fixed parameters
#' @param mc.maxpts number of Monte-Carlo samples
#'
#' @return the Monte-Carlo matrix, first column is always the 'sample' column
#' @engine
sampleMC_norm <- function(tdmorefit, regimen=NULL, fix=tdmorefit$fix, mc.maxpts=100) {
  regimen <- regimen %||% tdmorefit$regimen
  maxOcc <- getMaxOccasion(regimen)
  start <- processParameters(coef(tdmorefit), tdmorefit$tdmore, regimen)

  # prepare the variance-covariance matrix
  varcov <- expandOmega(tdmorefit$tdmore, maxOcc) #non-IOV terms are always first
  n <- nrow( vcov(tdmorefit) ) #number of occ available in tdmorefit
  n <- min( maxOcc, n) #if maxOcc is less, use that one
  varcov[1:n, 1:n] <- vcov(tdmorefit)[1:n, 1:n]

  parNames <- names(start)
  fixIndexes <- getFixedParametersIndexes(parNames = parNames, fix = fix)

  if(is.null(fixIndexes)) {
    mc <- mnormt::rmnorm(mc.maxpts, mean=start, varcov=varcov)
    if(mc.maxpts == 1) mc <- matrix(mc, nrow=1)
    mc <- as.data.frame(mc)
    colnames(mc) <- parNames
  } else {
    mc <- mnormt::rmnorm(mc.maxpts, mean=start[-fixIndexes], varcov=varcov[-fixIndexes, -fixIndexes])
    if(mc.maxpts == 1) mc <- matrix(mc, nrow=1)
    mc <- as.data.frame(mc)

    if(length(fix) > 1) {
      mc_fix <- t(replicate(n=nrow(mc), expr=fix))
    } else {
      mc_fix <- replicate(n=nrow(mc), expr=fix)
    }
    mc <- cbind(mc_fix, mc)
    colnames(mc) <- c(names(fix), parNames[-fixIndexes])
  }
  mc <- cbind( sample=1:mc.maxpts, mc )
  row.names(mc) <- NULL
  list(mc=mc)
}

#' Generate parameter samples using the Metropolis MCMC method.
#'
#' @param tdmorefit a tdmorefit object
#' @param fix named numeric vector with the fixed parameters
#' @param mc.maxpts number of Monte-Carlo samples
#' @param mc.batch number of points to run per batch
#' @param tune tuning parameter to reduce `omega`
#' @inheritParams MCMCpack::MCMCmetrop1R
#' @param .progress Allows to specify a plyr-like progress object
#' A plyr progress object is a list with 3 function definitions: `init(N)`, `step()` and `term()`.
#' This can also be specified as a boolean. TRUE uses the default dplyr progress_estimated.
#'
#' @return the Monte-Carlo matrix, first column is always the 'sample' column
#' @engine
sampleMC_metrop <- function(tdmorefit, regimen=NULL, fix=tdmorefit$fix, mc.maxpts=100, mc.batch=max(10, floor(mc.maxpts/2)), verbose=0, tune=1, .progress=interactive()) {
  if(!is.null(mc.maxpts)) {
    p <- to_progress(.progress)
    p$initialize(total=mc.maxpts)

    N <- 0
    sampled <- 0
    accepted <- 0
    result <- list()
    while(N < mc.maxpts) { #sample until we have enough
      out <- sampleMC_metrop(tdmorefit, fix=fix, mc.maxpts=NULL, mc.batch=mc.batch, verbose=verbose, tune=tune)
      # https://projecteuclid.org/euclid.aoap/1034625254
      # Tune the proposal variance so that the average acceptance rate is roughly 1/4.
      if(out$acceptance < 0.1) tune = tune / 2 #jump around less
      else if (out$acceptance > 0.5) tune = tune * 2 #jump around more

      N <- N + nrow(out$mc)
      p$tick(nrow(out$mc))
      result <- c(result, list(out) )
    }

    mc <- lapply(seq_along(result), function(i) {
      mc <- result[[i]]$mc
      if(nrow(mc) == 0) return(NULL)
      cbind(chain=i, result[[i]]$mc)
    })
    mc <- purrr::reduce(mc, rbind) #required instead of bind_rows, to maintain duplicate columns in case of IOV
    mc <- mc[seq(length.out=mc.maxpts), ] #only keep the first X rows, not all rows!
    ### right now, we have 'chain' and 'sample' columns
    ## avoid any dplyr verbs, because they cannot handle repeated column names
    parNames <- colnames(mc)[c(-1, -2)]
    mc <- cbind(sample=seq_len(nrow(mc)), chain=mc$chain, chain.sample=mc$sample, mc[, c(-1, -2), drop=FALSE]) # add chain.sample column
    colnames(mc)[4:ncol(mc)] <- parNames
    ## sample, chain, chain.sample, rest

    sampled=sum( purrr::map_dbl(result, ~.x$sampled) )
    accepted=sum( purrr::map_dbl(result, ~.x$accepted) )
    acceptance=accepted / sampled
    return( list(
      mc=mc,
      sampled=sampled,
      accepted=accepted,
      acceptance=acceptance
    ))
  }
  regimen <- regimen %||% tdmorefit$regimen
  maxOcc <- getMaxOccasion(regimen)
  start <- processParameters(coef(tdmorefit), tdmorefit$tdmore, regimen)

  ## Ok, thinking time
  ## If you use omega, it will generate proposal distributions that are only based on the initial population,
  ## BUT are possibly very far from the posthoc distribution.
  ## If this is the case, the MH-sampler will not accept any proposed parameters, since all generated proposals are
  ## unlikely in the posthoc distribution. Acceptance rate will be 0.
  ##
  ## In essence, we are jumping around a large parameter space with only a very small space
  ## that yields good options!
  ##
  ## An alternative is to jump around using the actual variance-covariance matrix from the optimization step.
  ## However, the reason we are in MH-sampling is exactly because this variance-covariance matrix is difficult to establish!
  ##
  ## A pragmatic solution is to use a smaller omega to jump around, increasing the acceptance rate. Yes,
  ## this is a biased estimate of uncertainty, but it is better than nothing...
  omega <- tdmorefit$tdmore$omega
  omega <- expandOmega(tdmorefit$tdmore, maxOcc)

  omegaChol <- Matrix::chol(omega)

  tdmore <- tdmorefit$tdmore
  cTdmore <- tdmore
  cTdmore$cache <- model_prepare(tdmore$model, times=tdmorefit$observed$TIME, regimen=tdmorefit$regimen, parameters=start, covariates=tdmorefit$covariates, iov=tdmore$iov, extraArguments=tdmore$extraArguments)

  fun <- function(par) {
    names(par) <- names(start)
    ll(par=par, omega=omegaChol, isChol=TRUE, fix=fix, tdmore=cTdmore, observed=tdmorefit$observed, regimen=tdmorefit$regimen, covariates=tdmorefit$covariates)
  }
  verboseOutput <- utils::capture.output({
    out <- MCMCpack::MCMCmetrop1R(fun, theta.init=start, seed=floor(stats::runif(1)*1E6), burnin=0, mcmc=mc.batch, V=omega*tune, verbose=verbose)
  })
  if(verbose > 0) cat(verboseOutput)

  #rejectionRate <- coda::rejectionRate(out)
  #effRate <- coda::effectiveSize(out) / mc.batch

  rejected <- c(FALSE, #first row is always OK!
                out[-nrow(out),1, drop=TRUE] == out[-1,1,drop=TRUE]) #just test on column 1; sufficient

  mc <- as.data.frame(out)
  colnames(mc) <- names(start)
  mc <- mc[ !rejected, , drop=F]
  mc <- cbind( sample=seq_len(nrow(mc)), mc )

  list(mc=mc, sampled=nrow(out), accepted=sum(!rejected), acceptance=sum(!rejected) / nrow(out) )
}

#' Predict new data using a model fit
#'
#' @param object A tdmorefit object
#' @param newdata
#' A data.frame with new data and the columns to predict,
#' or a numeric vector to specify times, and predict all model output
#' or NULL to interpolate between 0 and the maximum known times
#' @param regimen Treatment regimen
#' @param parameters named numeric vector of fixed parameters
#' @param covariates the model covariates, named vector, or data.frame with column 'TIME', and at least TIME 0
#' @param se.fit TRUE to provide a confidence interval on the prediction, adding columns xxx.median, xxx.upper and xxx.lower
#' FALSE to show the model prediction (IPRED)
#' @param level The confidence interval, or NA to return all mc.maxpts results
#' @param mc.maxpts Maximum number of points to sample in Monte Carlo simulation
#' @param ... ignored
#' @param .progress Allows to specify a plyr-like progress object
#' A plyr progress object is a list with 3 function definitions: `init(N)`, `step()` and `term()`.
#' This can also be specified as a boolean. TRUE uses the default dplyr progress_estimated.
#'
#' @return A data.frame
#'
#' @importFrom stats coef vcov median quantile
#' @export
predict.tdmorefit <- function(object, newdata=NULL, regimen=NULL, parameters=NULL, covariates=NULL, se.fit=FALSE, level=0.95, mc.maxpts=100, ..., .progress=interactive()) {
  tdmorefit <- object
  if(is.null(regimen)) regimen=tdmorefit$regimen
  if(is.null(newdata)) newdata <- model.frame(tdmorefit)

  par <- processParameters(parameters, tdmorefit$tdmore, regimen, defaultValues=coef(tdmorefit))
  if(is.null(covariates)) covariates <- tdmorefit$covariates

  ipred <- stats::predict(object=tdmorefit$tdmore, newdata=newdata, regimen=regimen, parameters=par, covariates=covariates)
  if(mc.maxpts > 0 && se.fit) {
    p <- to_progress(.progress)
    p$initialize(total=mc.maxpts)

    out <- sampleMC(tdmorefit, regimen=regimen, fix = parameters, mc.maxpts = mc.maxpts)
    mc <- out$mc
    ## remove 'chain' and 'chain.sample'
    if(colnames(mc)[2] == "chain" && colnames(mc)[3] == "chain.sample") {
      ## remove these values
      mc['chain'] <- NULL
      mc['chain.sample'] <- NULL
      #mc <- mc[ , c(-2, -3), drop=FALSE] ## This renames the columns!!
    }

    uniqueColnames <- make.unique(colnames(mc)) # needed for dplyr to have unique colnames

    # Prepare the tdmore cache
    cTdmore <- tdmorefit$tdmore
    cTdmore$cache <- model_prepare(cTdmore$model, times=ipred$TIME, regimen=regimen, parameters=par, covariates=covariates, iov=cTdmore$iov, extraArguments=cTdmore$extraArguments)
    fun <- function(i) {
      p$tick()
      row <- mc[i,,drop=TRUE] #make vector
      res <- unlist(row[-1]) # Remove 'sample'
      names(res) <- names(par)
      pred <- stats::predict(object=cTdmore, newdata=newdata, regimen=regimen, parameters=res, covariates=covariates)
      names(row) <- uniqueColnames
      if(nrow(pred)==0) row <- lapply(row, function(x){rep(x, length.out=0)}) #just bind the names, not the values
      resArray <- cbind(row, pred)
      resArray
    }
    fittedMC <- lapply(mc$sample, fun)
    fittedMC <- purrr::reduce(fittedMC, rbind) #required instead of bind_rows, to maintain duplicate columns in case of IOV
    colnames(fittedMC)[seq_len(length(uniqueColnames))] <- colnames(mc)

    if(is.na(level)) { #user requested full dataframe without summary
      return(fittedMC)
    }
    oNames <- getPredictOutputNames(newdata, colnames(fittedMC), names(coef(tdmorefit)))
    summary <- summariseFittedMC(fittedMC, ipred, level, oNames)
    return(summary)
  } else if (mc.maxpts==0 && se.fit) {
    oNames <- getPredictOutputNames(newdata, colnames(ipred), names(coef(tdmorefit)))
    summary <- summariseFittedMC(ipred, ipred, level, oNames)
    return(summary)
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
#' @engine
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
#' @engine
summariseFittedMC <- function(fittedMC, ipred, level, oNames) {
  if(length(oNames)==0) {
    return( fittedMC %>% dplyr::distinct(.data$TIME) )
  }
  a <- (1-level)/2
  #dplyr gets confused with multiple 'ECL' columns (in case of IOV)
  fittedMC <- fittedMC[,c("TIME", oNames)]

  lower <- function(x) {quantile(x, probs=a)}
  upper <- function(x) {quantile(x, probs=1-a)}
  retValue <- fittedMC %>%
    dplyr::group_by(.data$TIME) %>%
    dplyr::summarize_at(dplyr::vars(!!oNames), list(median, lower, upper)) %>%
    dplyr::ungroup()
  ## Assign the right names to these columns
  names <- expand.grid(oNames, c("median", "lower", "upper"))
  cNames <- paste0(names$Var1, ".", names$Var2)
  colnames(retValue)[seq_along(cNames)+1] <- cNames

  retValue[, oNames] <- ipred[, oNames]

  return(retValue)
}

#' Calculate the log-likelihood of the predicted values.
#' If fixed values were used to obtain the tdmorefit, they are considered estimated and are included in the population log-likelihood.
#'
#' @param object A tdmorefit object
#' @param type log-lokelihood function type, 3 possible values: 'pop', 'pred' or 'll' (= pop + pred)
#' @param par the parameter set to evaluate in log likelihood
#' @param ... ignored
#'
#' @return A numeric value
#'
#' @importFrom stats formula model.frame logLik
#' @export
logLik.tdmorefit <- function(object, type=c('ll', 'pop', 'pred'), par=coef(object), ...) {
  tdmore <- formula(object)
  observed <- model.frame(object)
  regimen <- object$regimen
  covariates <- object$covariates
  fun <- getLikelihoodFun(type)
  omega <- expandOmega(tdmore, getMaxOccasion(regimen))

  fun(par=par,
      omega=omega,
      fix=list(indexes=c(), values=c()),
      tdmore=tdmore,
      observed=observed,
      regimen=regimen,
      covariates=covariates)
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
#' @param ... ignored
#'
#' @inheritParams model.frame.tdmore
#'
#' @return
#' A data.frame, similar to the one used to estimate this `tdmorefit` object.
#' If `se` was specified, then a column xx.upper and xx.lower with the
#' lower and upper confidence interval (based on the residual error model) is added.
#'
#' @export
model.frame.tdmorefit <- function(formula, data=NULL, se=FALSE, level=0.95, onlyOutput=FALSE, ...) {
  if(is.null(data)) data <- formula$observed
  result <- model.frame.tdmore(formula=formula$tdmore, data=data, se=se, level=level, onlyOutput=onlyOutput, ...)
  result
}

#' Get the right likelihood function for the specified type.
#'
#' @param type likelihood type, ll, pop or pred
#'
#' @return the correct tdmore function
#' @noRd
getLikelihoodFun <- function(type) {
  type <- match.arg(type, c('ll', 'pop', 'pred'))
  switch(type,
         ll=ll,
         pop=pop_ll,
         pred=pred_ll
  )
}

#' Check if the given argument is of class `tdmorefit`
#' @param x An object to test
#' @keywords internal
#' @engine
is.tdmorefit <- function(a) {inherits(a, "tdmorefit")}

#' Expand OMEGA matrix by duplicating variance-covariance information for IOV terms.
#' Note that the column and row names of the returned matrix are strictly identical to the ones returned by getParameterNames().
#'
#' @param tdmore the tdmore object
#' @param occasions how many occasions
#' @return the omega matrix with the duplicated IOV terms at the end
#' @noRd
expandOmega <- function(tdmore, occasions) {
  iov <- tdmore$iov
  omega <- tdmore$omega

  if(is.null(iov)) {
    return(omega)
  }

  iov <- colnames(omega) %in% iov

  elements <- list(omega[!iov, !iov], omega[iov, iov] )
  result <- Matrix::.bdiag(
    elements[c(1, rep(2, occasions) ) ]
  ) %>% as.matrix()

  myNames <- rownames(tdmore$omega)
  myNames <- c(myNames[!iov], rep(myNames[iov], occasions) )
  dimnames(result) <- list(myNames, myNames)

  result
}

#' Get the residuals for a specific tdmorefit
#'
#' @param data optional data.frame
#' If specified, this adds the original residuals back onto this data.frame
#' If points at a time without corresponding IRES are present, they are not modified.
#'
#' @rdname residuals
#'
#' @examples
#' observed <- data.frame(TIME=5, CONC=0.060)
#' regimen=data.frame(TIME=0, AMT=10)
#' fit <- estimate(getModel("default"),
#'     regimen=regimen,
#'     observed=observed,
#'     covariates=c(WT=70)
#'     )
#' wres <- residuals(fit, weighted=TRUE)
#'
#' all.equal(
#'   residuals(fit, predict(fit), weighted=TRUE),
#'   observed
#' )
#'
#' regimen <- fit$regimen
#' regimen$AMT <- regimen$AMT*2
#' predictionForDoubleDose <- residuals(fit, predict(fit, regimen=regimen), weighted=TRUE)
#'
#' @importFrom stats residuals
#' @noRd
#' @export
residuals.tdmorefit <- function(object, data, weighted=FALSE, ...) {
  res <- residuals(object$tdmore, stats::predict(object), model.frame(object), weighted=weighted, ...)
  if(missing(data)) return(res)

  observed <- model.frame(object)
  if(any(duplicated(observed$TIME)) || any(duplicated(data$TIME)) ) stop()
  i <- data$TIME %in% observed$TIME #matching rows from data
  j <- observed$TIME %in% data$TIME #matching rows from observed

  result <- residuals.tdmore(object$tdmore, predicted=data[i,], observed=res[j,], weighted=TRUE, inverse=TRUE)
  data[i, names(result)] <- result[i, ]
  return(data)
}

#' Create a typical value subject
#'
#' @param x a tdmorefit object, or tdmore object
#' @param ... extra parameters passed to the internal engine
#' @export
as.population <- function(x, ...) {
  estimate(x, observed=model.frame(x, data=numeric()), ...)
  #estimate(x, observed=NA, control=list(trace=0))
}

#' Create a tiible of N random tdmorefit instances
#' @param x a tdmorefit or tdmore object
#' @param N the number of random subjects to create
#'
#' @details This method samples ETA parameters from the variance-covariance
#' matrix of a given tdmorefit object. If a tdmore object is provided,
#' the population is used.
#'
#' The regimen needs to be provided to know the number of a
#'
#' @export
as.sample <- function(x, N=100) {
  if(is.tdmore(x)) x <- as.population(x)
  if(!is.tdmorefit(x)) stop("Please provide a tdmorefit object")
  par <- x$res
  out <- sampleMC(x, mc.maxpts=N)
  mc <- out$mc
  colSelect <- colnames(mc) %in% names(par)
  varcov <- diag(1e-16, length(par), length(par))
  result <- purrr::map(seq_len(nrow(mc)), function(i){
    candidate <- x
    row <- mc[i,,drop=TRUE] #make vector
    res <- unlist(row[colSelect]) # Remove 'sample'
    names(res) <- names(par)
    candidate$res <- res
    candidate$varcov <- varcov

    candidate
  })
  tibble(ID = seq_len(N),
         fit=result
  )
}
