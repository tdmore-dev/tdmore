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
#' @param model
#' @param observed
#' @param ...
#' Extra parameters to pass to nlm
#'
#' @return
#' @export
#'
#' @examples
estimate <- function(tdmore, observed=data.frame(), regimen, ...) {
  assert_that(class(tdmore) == "tdmore")
  pointEstimate <- nlm(f=function(...) {-2*ll(...)},
                       p=rep(0, length(tdmore$parameters)),
                       tdmore,
                       observed,
                       regimen,
                       print.level=2,
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
      varcov=varcov,
      model=model
    ), class=c("tdmorefit"))
}


# Good generics
#    profile         residuals
#terms
summary.tdmorefit <- function(tdmorefit) {
  cat("Fit of model to ", nrow(tdmorefit$observed), " points of data\n")
  cat("OFV = ", tdmorefit$ofv, "\n")
  cat("Point estimate: ", tdmorefit$res, "\n")
  cat("Variance/covariance matrix: ", tdmorefit$varcov, "\n")
}

coef.tdmorefit <- function(tdmorefit) {tdmorefit$res}
vcov.tdmorefit <- function(tdmorefit) {tdmorefit$varcov}
confint.tdmorefit <- confint.default

fitted.tdmorefit <- function(tdmorefit) {
  predict.tdmore(tdmorefit$tdmore, tdmorefit$observed, tdmorefit$regimen, tdmorefit$res)
}

predict.tdmorefit <- function(tdmorefit, newdata, se.fit=FALSE, level=0.95) {
  ## TODO: error checking of newdata?
  assert_that("TIME" %in% colnames(newdata))
  predict.tdmore(tdmorefit$tdmore, newdata, tdmorefit$regimen, tdmorefit$res, se=se.fit)
}

logLik.tdmorefit <- function(tdmorefit) {tdmorefit$logLik}

profile.tdmorefit <- function(tdmorefit, which = 1:npar, maxpts = 100, limits=c(-1.96, 1.96)) {
  model <- tdmorefit$tdmore
  x <- seq(limits[1], limits[2], length.out=maxpts)
  x <- expand.grid(rep(list(x), length(model$parameters)))
  colnames(x) <- model$parameters
  result <- plyr::adply(x, 1, function(estimate) {
    c(logLik=ll(as.numeric(estimate), model, tdmorefit$observed, tdmorefit$regimen))
  }, .progress="text")
  result
}




## TODO: add an argument 'interval' to predict and 'maxpts'
predictMC.EstimationResult <- function(estimationResult, times, N=100) {
  mc <- mnormt::rmnorm(N, mean=estimationResult$res, varcov=estimationResult$varcov) %>% as.data.frame()
  model <- estimationResult$model
  colnames(mc) <- model$parameters
  mc$subject <- 1:N
  fittedMC <- plyr::ddply(mc, 1, function(res) {
    res <- res[names(res) != "subject"]
    pred <- predict.TDM.Model(model, estimates=res, times=times) %>% as.data.frame
    pred$subject <- res$subject
    pred
  }, .progress="text")
  fittedMC
}
