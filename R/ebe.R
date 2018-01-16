# The functions here are used to calculate an emperical bayesian estimate
pop_ll <- function(estimate) {
  sum( dnorm(estimate, log=TRUE ) )
}

pred_ll <- function(estimate, tdmore, observed, regimen) {
  if(nrow(observed) == 0) return(0)
  pred <- predict.tdmore(tdmore, observed, regimen, estimate)
  res <- residuals.tdmore(tdmore, observed, pred, log=TRUE)
  sum(res)
}

ll <- function(estimate, tdmore, observed, regimen) {
  names(estimate) <- tdmore$parameters
  res <- pop_ll(estimate) + pred_ll(estimate, tdmore, observed, regimen)
  res
}


#' Calculate the Emperical Bayesian Estimated parameters that predict
#' a given dataset using the model
#'
#' @param model a tdmore object
#' @param observed data frame with at least a TIME column, and all observed data. The observed data will be compared to the model predictions.
#' @param regimen data frame
#' @param ... Extra parameters to pass to nlm
#'
#' @return
#' @export
#'
#' @examples
estimate <- function(tdmore, observed=data.frame(), regimen, print.level=0, ...) {
  assert_that(class(tdmore) == "tdmore")
  pointEstimate <- nlm(f=function(...) {-2*ll(...)},
                       p=rep(0, length(tdmore$parameters)),
                       tdmore,
                       observed,
                       regimen,
                       print.level=print.level,
                       hessian=TRUE,
                       ...)
  res <- pointEstimate$estimate
  names(res) <- tdmore$parameters

  #OFIM <- pointEstimate$hessian ##!!! Wrong hessian!?
  # Observed fisher information matrix is -1*hessian(ll)
  OFIM <- numDeriv::hessian(ll, res, tdmore=tdmore, observed=observed, regimen=regimen)
  varcov <- solve(OFIM) #inverse of OFIM is an estimator of the asymptotic covariance matrix
  if(all(diag(varcov) <= 0)) varcov <- -1 * varcov
  dimnames(varcov) = list(tdmore$parameters, tdmore$parameters)
  ofv <- pointEstimate$minimum

  structure(
    list(
      tdmore=tdmore,
      observed=observed,
      regimen=regimen,
      ofv=ofv,
      logLik=ofv / -2,
      res=res,
      varcov=varcov
    ), class=c("tdmorefit"))
}

#' Title
#'
#' @param tdmorefit
#'
#' @return
#' @export
#' @S3method
#'
#' @examples
summary.tdmorefit <- function(tdmorefit) {
  cat("Fit of model to ", nrow(tdmorefit$observed), " points of data\n")
  cat("LogLikelihood = ", logLik(tdmorefit), "\n")
  cat("Point estimate: ", coef(tdmorefit), "\n")
  cat("Variance/covariance matrix: ", vcov(tdmorefit), "\n")
}

#' Title
#'
#' @param tdmorefit
#'
#' @return
#' @export
#' @S3method
#'
#' @examples
coef.tdmorefit <- function(tdmorefit) {tdmorefit$res}
#' Title
#'
#' @param tdmorefit
#'
#' @return
#' @S3method
#' @export
#'
#' @examples
vcov.tdmorefit <- function(tdmorefit) {tdmorefit$varcov}

#' @S3method
confint.tdmorefit <- function(tdmorefit) {
  confint.default(tdmorefit)
}

#' @S3method
fitted.tdmorefit <- function(tdmorefit) {
  predict.tdmore(tdmorefit$tdmore, tdmorefit$observed, tdmorefit$regimen, tdmorefit$res)
}


#' Title
#'
#' @param tdmorefit
#' @param newdata
#' @param se.fit
#' @param level
#' @param maxpts
#'
#' @return
#' @S3method
#' @export
#'
#' @examples
predict.tdmorefit <- function(tdmorefit, newdata, se.fit=FALSE, level=0.95, maxpts=100, .progress="none") {
  assert_that("TIME" %in% colnames(newdata))

  ipred <- predict.tdmore(tdmorefit$tdmore, newdata, tdmorefit$regimen, tdmorefit$res)
  if(se.fit) {
    oNames <- names(newdata)
    oNames <- oNames[oNames != "TIME"]
    pNames <- names(coef(tdmorefit))
    mc <- as.data.frame( mnormt::rmnorm(maxpts, mean=coef(tdmorefit), varcov=vcov(tdmorefit)) )
    mc$subject <- 1:maxpts
    fittedMC <- plyr::ddply(mc, 1, function(row) {
      res <- row[pNames]
      pred <- predict.tdmore(tdmorefit$tdmore, newdata, tdmorefit$regimen, unlist(res))
      pred$subject <- row$subject
      pred
    }, .progress=.progress)
    a <- (1-level)/2
    plyr::ddply(fittedMC, plyr::.(TIME), function(x) {
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

#' @S3method
logLik.tdmorefit <- function(tdmorefit) {tdmorefit$logLik}

#' @S3method
profile.tdmorefit <- function(tdmorefit, which = 1:npar, maxpts = 100, limits=c(-1.96, 1.96), .progress="none") {
  model <- tdmorefit$tdmore
  x <- seq(limits[1], limits[2], length.out=maxpts)
  x <- expand.grid(rep(list(x), length(model$parameters)))
  colnames(x) <- model$parameters
  result <- plyr::adply(x, 1, function(estimate) {
    c(logLik=ll(as.numeric(estimate), model, tdmorefit$observed, tdmorefit$regimen))
  }, .progress=.progress)
  result
}

formula.tdmorefit <- function(tdmorefit) {tdmorefit$model}

model.frame.tdmorefit <- function(tdmorefit, se=FALSE, level=0.95) {
  model.frame.tdmore(tdmorefit$tdmore, tdmorefit$observed, se=se, level=level)
}
