#' MPC predict method. This simply generates a cumulative sum for all
#'
#' @inheritParams predict.tdmore
#' @return a data.frame with all observed values at the given time points
#' @engine
predict.tdmore_mpc <- function(object, newdata, regimen=NULL, parameters=NULL, covariates=NULL, se=FALSE, level=0.95, ...) {
  for(par in object$iov) {
    i <- which(names(parameters) == par)
    parameters[i] <- cumsum(parameters[i])
  }

  NextMethod()
}

#' MPC estimate method.
#'
#' @inheritParams estimate
#' @param .progress progress bar definition
#' @return a tdmorefit
#'
#' @importFrom dplyr filter
#' @engine
estimate_tdmore_mpc <- function(object, observed=NULL, regimen=NULL, covariates=NULL, par=NULL, fix=NULL, method="L-BFGS-B", se.fit=TRUE, lower=NULL, upper=NULL, multistart=F, control=list(), ..., .progress=NULL) {
  p <- to_progress(.progress)
  objectNoMpc <- object
  class(objectNoMpc) <- setdiff(class(object), "tdmore_mpc")

  if(!"OCC" %in% colnames(regimen)) {
    warning("OCC column missing, adding...")
    regimen$OCC <- seq_len(nrow(regimen))
  }
  stopifnot(is.null(fix)) #MPC assumes FIX is empty
  if(nrow(regimen)==0 || is.null(observed) || nrow(observed)==0) {
    fit <- estimate(
      objectNoMpc,
      observed=observed, regimen=regimen, covariates=covariates, par=par, fix=fix,
      method=method, se.fit=se.fit, lower=lower, upper=upper, multistart=multistart,
      control=control, ..., .progress=.progress
    )
    fit$tdmore <- object
    return(fit)
  }

  # advance the fit per occasion
  occasions <- tibble(
    OCC=unique(regimen$OCC)
    )
  occasions$from <- regimen$TIME[ match(occasions$OCC, regimen$OCC) ] #shows first match
  occasions$from[1] <- 0
  occasions$to <- c(occasions$from[-1], Inf)

  fix <- c()
  estim <- function(from, to) {
    thisRegimen <- filter(regimen, .data$TIME < to) #strictly less !
    #only observations in the current occasion (and before, since those parameters cannot move anymore anyway)
    #thisObserved <- filter(observed, .data$TIME <= to) #we count ON the occasion as well
    thisObserved <- filter(observed, .data$TIME > from & .data$TIME <= to) #we count ON the occasion as well
    fit <- estimate(
      objectNoMpc, observed=thisObserved, regimen=thisRegimen, covariates=covariates, par=par, fix=fix, method=method, se.fit=se.fit,
      lower=lower, upper=upper, multistart=multistart, control=control, ...
    )
    fix <<- coef(fit)
    fit
  }
  p$initialize(total=length(occasions$to))
  fits <- lapply(1:nrow(occasions), function(i){
    p$tick()
    estim(occasions$from[i], occasions$to[i])
  })

  fit <- fits[[ length(fits) ]] #last fit
  fit$res <- coef(fit)
  fit$fix <- c()
  varcovsMissing <- vapply(fits, function(x) is.null(x$varcov), FUN.VALUE=TRUE)
  if(any(varcovsMissing)) {
    fit$varcov <- NULL
  } else {
    varcovs <- lapply(fits, function(x) {
      i <- seq_len(nrow(x$varcov))
      if(length(x$fix) > 0) i <- utils::tail(i, n=-1*length(x$fix)) #only the actual varcov, not the FIX part
      x$varcov[i, i]
    })
    fit$varcov <- Matrix::.bdiag(varcovs) %>% as.matrix()
  }
  fit$observed <- observed
  fit$tdmore <- object

  return(fit)
}

#' MPC is a generic function to make a tdmore model compatible with MPC.
#'
#' @param x a model object (type not defined)
#' @param parameters select which parameters should be mpc parameters. MPC parameters are cumulatively summed across occasions.
#' @param ... extra arguments
#' @return an object of class tdmore_mpc
#' @export
mpc <- function(x, parameters, ...) {
  UseMethod("mpc")
}

#' @export
mpc.tdmore <- function(x, parameters, ...) {
  class(x) <- append("tdmore_mpc", "tdmore")
  if( !setequal(x$parameters, x$iov) ) stop("For MPC, all parameters should be IOV")
  return(x)
}

#'
#' Test if an object is a `tdmore_mpc` model
#'
#' @param x object to test
#' @return TRUE if x inherits from `tdmore_mpc`
#' @noRd
is.mpc <- function(x) {
  inherits(x, "tdmore_mpc")
}

#' @export
print.tdmore_mpc <- function(x, ...) {
  cat("Mpc model\n")
  NextMethod()
}
