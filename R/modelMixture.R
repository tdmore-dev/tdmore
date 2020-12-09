#' Create a TDMore mixture model.
#'
#' @param ... 2 or more tdmore models that describe different subpopulations
#' @param probs 'a priori' probabilities for belonging to the different subpopulations, numeric vector
#'
#' @return a tdmore_mixture object
#' @export
tdmore_mixture <- function(..., probs) {
  models <- list(...)
  for (model in models) {
    if(!is.tdmore(model)) stop("Only tdmore models can be added to the tdmore mixture")
  }
  if(! length(models) >= 2) stop("You should provide at least two tdmore models to create a tdmore mixture")
  if(!is.numeric(probs)) stop("probs is not numeric")
  if(length(models) != length(probs)) stop("There must be as many models as probabilities defined in the numeric vector 'probs'")
  if(!isTRUE(all.equal(target = 1, current = sum(probs), tolerance = 1e-6))) stop("Sum of probabilities must be 1")
  tdmoreMixture <- structure(list(
    models=models,
    probs=probs
  ), class="tdmore_mixture")

  return(tdmoreMixture)
}

#' Create a TDMore mixture model, based on several options for a discrete covariate.
#'
#' @param model a tdmore model with a discrete covariate
#' @param probs a numeric vector with the a priori probabilities for the covariate values
#' @param covariates a list with covariate values corresponding to the right names
#'
#' @return a tdmore_mixture object
#'
#' @export
#' @importFrom purrr map
#' @examples
#' \dontrun{
#' tdmore_mixture_covariates(tacrolimus_storset,
#'           probs=c(0.5, 0.5),
#'           covariates=list(c(CYP3A5=1), c(CYP3A5=0))
#' )
#' }
tdmore_mixture_covariates <- function(model, probs, covariates) {
  N <- length(probs)
  stopifnot(length(covariates) == N)
  models <- purrr::map(covariates, ~ curry_covariate(model, .x) )
  do.call(tdmore_mixture, c(models, list(probs=probs))  )
}


#' Fill in a specific covariate in a model, transforming it into a more specific model.
#'
#' @param model tdmore model with a discrete covariate
#' @param curriedCovariates the curried covariates
#'
#' @return the model with the curried covariates
#' @noRd
curry_covariate <- function(model, curriedCovariates) {
  stopifnot( length(names(curriedCovariates)) == length(curriedCovariates) ) #all covariates should be named
  stopifnot( all(names(curriedCovariates) %in% model$covariates) )
  stopifnot(  all(is.numeric(curriedCovariates))  )
  class(model) <- c("tdmoreCurried", class(model))
  model$curriedCovariates <- curriedCovariates
  model
}

#' Predict a tdmoreCurried object.
#'
#' @param object a tdmoreCurried object
#' @param ... all the extra arguments to pass to the next predict method
#'
#' @return a dataframe with the predictions
#' @engine
predict.tdmoreCurried <- function(object, ...) {
  args <- list(...)

  ## Covariates could be data.frame or numeric vector
  if(is.null(args$covariates)) {
    covariates <- object$curriedCovariates
  } else {
    if(is.data.frame(args$covariates)) {
      covariates <- cbind( args$covariates, object$curriedCovariates )
    } else {
      covariates <- c(object$curriedCovariates, args$covariates)
    }
  }
  class(object) <- setdiff(class(object), "tdmoreCurried" ) # Remove tdmoreCurried class
  NextMethod(object=object, covariates=covariates, ...)
}

#' @export
print.tdmoreCurried <- function(x, ...) {
  cat("Pre-filled covariates: \n")
  print(x$curriedCovariates)
  NextMethod()
}

#' Predict from a TDMore mixture model. The model with the highest probability is used for the predictions.
#'
#' @inheritParams predict.tdmore
#' @param ... extra arguments for the call to the model
#'
#' @return A data.frame with all observed values at the given time points
#' @engine
predict.tdmore_mixture <- function(object, newdata, regimen=NULL, parameters=NULL, covariates=NULL, se=FALSE, level=0.95, ...) {
  tdmoreMixture <- object
  chosenModel <- tdmoreMixture$models[[which.max(tdmoreMixture$probs)]]
  return(chosenModel %>% stats::predict(newdata, regimen, parameters, covariates, se, level, ...))
}

#' Calculate the Empirical Bayesian Estimate parameters that predict
#' a given dataset using a mixture model.
#'
#' @inheritParams estimate
#' @param ... extra parameters to pass to the optim function
#'
#' @return A tdmorefit_mixture object
#' @importFrom stats optim
#' @engine
estimateMixtureModel <- function(object, observed=NULL, regimen, covariates=NULL, par=NULL, fix=NULL, method="L-BFGS-B", lower=NULL, upper=NULL, ...) {
  mixture <- object
  models <- mixture$models
  probs <- mixture$probs
  fits <- lapply(models, FUN=function(model) {
    estimate(model, observed=observed, regimen=regimen, covariates=covariates, par=par, fix=fix, method=method, lower=lower, upper=upper, ...)
  })
  mixtureProbs <-
    data.frame(lik = as.numeric(lapply(fits, FUN = function(fit) {exp(fit$logLik)})),
               probs = probs)
  mixtureProbs$IPkNumerator <- mixtureProbs$lik * mixtureProbs$probs
  mixtureProbs$IPk <- mixtureProbs$IPkNumerator/sum(mixtureProbs$IPkNumerator)

  retValue <- structure(list(
    mixture = mixture,
    fits = fits,
    fits_prob = mixtureProbs,
    winner = which(mixtureProbs$IPk == max(mixtureProbs$IPk))[1] # [1] in case of several winners
  ), class = c("tdmorefit_mixture", "tdmorefit"))

  return(retValue)
}

#' Predict from a tdmorefit mixture object.
#'
#' @inheritParams predict.tdmorefit
#' @param ... ignored
#'
#' @return a data.frame
#' @engine
predict.tdmorefit_mixture <- function(object, newdata=NULL, regimen=NULL, parameters=NULL, covariates=NULL, se.fit=FALSE, level=0.95, mc.maxpts=100, ...) {
  fits <- object$fits
  mixture <- object$mixture
  winner <- object$winner
  ipred <- stats::predict(fits[[winner]], newdata=newdata, regimen=regimen, parameters=parameters, covariates=covariates, se.fit=F)

  if(se.fit==TRUE) {
    fits_prob <- object$fits_prob
    mixnum <- sample(seq_along(fits_prob$IPk), size=mc.maxpts, replace=T, prob=fits_prob$IPk)
    fits_prob <- cbind(fit=seq_len(nrow(fits_prob)), fits_prob) # Fit column added with fit index
    fittedMC <- purrr::map_dfr(fits_prob$fit, function(i) {
          row <- fits_prob[i, , drop=FALSE]
              fitIndex <- row[["fit"]]
              samples <- sum(mixnum==fitIndex)
              fit <- fits[[fitIndex]]
             stats::predict(fit, newdata=newdata, regimen=regimen, parameters=parameters, covariates=covariates, se.fit=T, level=NA, mc.maxpts=samples)
           })

    if(is.na(level)) { #user requested full dataframe without summary
      return(fittedMC)
    }
    oNames <- getPredictOutputNames(newdata, colnames(fittedMC), names(coef(fits[[1]])))
    retValue <- summariseFittedMC(fittedMC, ipred, level, oNames)
    return(retValue)

  } else {
    ipred
  }
}

