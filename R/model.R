assert_that <- assertthat::assert_that
# Structural model: how to predict --------------------------------------------------------
#' Prototype predict function
#'
#' @param model The model itself.
#' @param observed dataframe with at least a column 'TIME' and other values. A prediction will be generated for each filled-in value.
#' @param regimen dataframe with column 'TIME' and adhering to standard NONMEM specifications otherwise (columns AMT, RATE, CMT)
#' @param parameters dataframe with column 'TIME' and a column for each covariate and parameter
#'
#' @return
#' A data.frame similar to the observed data frame, but with predicted values.
#'
model_predict <- function(model, observed, regimen, parameters) {
  UseMethod("model_predict")
}

# Main API ----------------------------------------------------------------
#' Create a new TDM model
#'
#' @param model
#' @param parameters character vector of parameter names, or NULL to auto-detect
#'
#' @return
#' a tdmore object
#' @export
#'
#' @examples
tdmore <- function(model, parameters=NULL, add=0, prop=0, exp=0, ...) {
  UseMethod("tdmore")
}

tdmore.default <- function(model, ...) {
  stop("Model of class ", class(model), " not supported.")
}

predict.tdmore <- function(tdmore, observed, regimen, parameters, se=FALSE, level=0.95) {
  predicted <- model_predict(tdmore$model, observed, regimen, parameters)
  model.frame.tdmore(tdmore, predicted, se=se, level=level)
}

residuals.tdmore <- function(tdmore, observed, predicted, log=TRUE) {
  oNames <- colnames(observed)
  oNames <- oNames[oNames != "TIME"]
  i <- oNames %in% colnames(predicted)
  i <- oNames[i]
  ipred <- predicted[,i]
  obs <- observed[,i]
  res <- tdmore$res_var
  if(res$exp != 0) {
    sd <- res$exp
    dnorm(log(ipred), log(obs), sd, log=log) #calculate residual in the log domain
  } else {
    sd <- sqrt(res$add**2 + (res$prop * ipred)**2)
    ifelse(sd == 0 & ipred==obs, #treat points with 0 res. error special
           0, #ignore point
           dnorm(ipred, obs, sd, log=log) # calculate residual
    )
  }
}

model.frame.tdmore <- function(tdmore, observed, se=FALSE, level=0.95) {
  if(!se) return(observed)
  oNames <- colnames(observed)
  oNames <- oNames[oNames != "TIME"]
  a <- (1 - level) / 2
  q <- qnorm(a)
  res <- tdmore$res_var
  for(n in oNames) {
    obs <- observed[, n]
    if(res$exp != 0) {
      observed[, paste0(n, ".upper")] <- obs * exp(res$exp * q)
      observed[, paste0(n, ".lower")] <- obs * exp(-res$exp * q)
    } else {
      observed[, paste0(n, ".upper")] <- obs * (1 + res$prop*q) + res$add*q
      observed[, paste0(n, ".lower")] <- obs * (1 - res$prop*q) - res$add*q
    }
  }

  observed
}

formula.tdmore <- function(tdmore) {tdmore$model}
