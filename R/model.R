#' Create a new TDM model
#'
#' @param predict
#' Function to execute the model prediction, corresponding to the definition
#' function(estimates, times)
#' The 'estimates' parameter should be, by default, a named numeric vector.
#' The 'times' parameter should be a numeric vector
#' This function should return a data frame with all predicted values
#' @param addSigma
#' @param propSigma
#'
#' @return
#' a TDM.Model object
#' @export
#'
#' @examples
Model <- function(predict, parameters, addSigma=0, propSigma=0) {
  structure(list(
    predict=predict,
    parameters=parameters,
    addSigma=addSigma,
    propSigma=propSigma
  ), class="TDM.Model")
}

predict.re.ci <- function(model, ...) {
  UseMethod("predict.re.ci")
}

# predict <- function(model, ...) {
#   UseMethod("predict")
# }

predict.TDM.Model <- function(model, times=NULL, estimate=NULL, oName="CONC", re=FALSE, conf.int=NULL) {
  predicted <- model$predict(times, estimate)
  if(re) {
    ipred <- predicted[,oName]
    sd <- sqrt(model$addSigma**2 + (model$propSigma*ipred)**2)
    ipredre <- rnorm(length(ipred), ipred, sd)
    ipredre[sd==0] <- ipred[sd==0]
    predicted[, paste0(oName, ".RE")] <- ipredre
  }
  if(!is.null(conf.int)) {
    ipred <- predicted[,oName]
    sd <- sqrt(model$addSigma**2 + (model$propSigma*ipred)**2)
    p <- (1 - conf.int)/2
    predicted[, paste0(oNames, ".CIupper")] <- ipred + sd * qnorm(p)
    predicted[, paste0(oNames, ".CIlower")] <- ipred - sd * qnorm(p)
  }
  predicted
}

residuals.TDM.Model <- function(model, estimate, observed, log=TRUE) {
  oNames <- colnames(observed)
  oNames <- oNames[oNames != "time"]
  i <- oNames %in% colnames(predicted)
  i <- oNames[i]
  ipred <- predicted[,i]
  obs <- observed[,i]
  sd <- sqrt(model$addSigma**2 + (model$propSigma*ipred)**2)
  ifelse(sd == 0 & ipred==obs, #treat points with 0 res. error special
         0, #ignore point
         dnorm(ipred, obs, sd, log=log) # calculate residual
  )
}
