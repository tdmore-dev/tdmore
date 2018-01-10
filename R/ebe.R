# The functions here are used to calculate an emperical bayesian estimate

pop_ll <- function(estimate) {
  sum( dnorm(estimate, log=TRUE ) )
}

pred_ll <- function(estimate, model, observed) {
  if(nrow(observed) == 0) return(0)
  pred <- predict.TDM.Model(model, estimate, observed$time)
  res <- residual.TDM.Model(model, pred, observed)
  sum(res)
}

ll <- function(estimate, model, observed) {
  names(estimate) <- model$parameters
  res <- pop_ll(estimate) + pred_ll(estimate, model, observed)
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
estimate.TDM.Model <- function(model, observed, ...) {
  pointEstimate <- nlm(f=function(...) {-1*ll(...)},
                       p=rep(0, length(model$parameters)),
                       model, observed,
                       print.level=2,
                       hessian=TRUE,
                       ...)
  res <- pointEstimate$estimate
  names(res) <- model$parameters
  #OFIM <- pointEstimate$hessian ##!!! Wrong hessian!?
  # Observed fisher information matrix is -1*hessian(ll)
  OFIM <- numDeriv::hessian(ll, res, model=model, observed=observed)
  varcov <- solve(OFIM) #inverse of OFIM is an estimator of the asymptotic covariance matrix
  if(all(diag(varcov) <= 0)) varcov <- -1 * varcov
  ofv <- pointEstimate$minimum

  structure(
    list(
      ofv=ofv,
      res=res,
      varcov=varcov,
      model=model
    ), class=c("EstimationResult"))
}

predict.EstimationResult <- function(estimationResult, times) {
  estimationResult$model$predict(
    estimationResult$res,
    times
  )
}

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

llSurface <- function(model, observed, limits=c(-1.96, 1.96), N=50) {
  x <- seq(limits[1], limits[2], length.out=N)
  x <- expand.grid(rep(list(x), length(model$parameters)))
  colnames(x) <- model$parameters
  plyr::adply(x, 1, function(estimate) {
    ll(as.numeric(estimate), model, observed)
  }, .progress="text")
}


ebefit <- function(x, ...) UseMethod("ebefit", x)
ebefit.default <- function(x, ...) {
  stop("Not possible to provide an ebefit for a model of class ", class(x))
}

