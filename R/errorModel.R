#' Instantiate a new error model.
#'
#' @param var the related variable name
#' @param add additive residual error, as stdev
#' @param prop proportional residual error, as stdev
#' @param type type of error: constant, proportional, combined1, combined2, exponential
#' or NULL to auto-detect. A combination of additive and proportional always assumes 'combined2'.
#'
#' @details Please note that convention differ between modeling software. While NONMEM would always code a
#' combined additive/proportional model as `W=SQRT( SIGMA(1)**2 + SIGMA(2)**2 * IPRED**2 )`,
#' Monolix uses a default of `CONC = Cc + (a + b*Cc)*e`. TDMore uses the NONMEM convention by default.
#'
#' @export
errorModel <- function(var="CONC", add=0, prop=0, type=NULL) {
  stopifnot(is.character(var))
  stopifnot(length(var)==1)
  if(add!=0 && prop==0) allowedType="constant"
  else if(prop!=0 && add==0) allowedType="proportional"
  else if (prop != 0 && add != 0) allowedType=c("combined2", "combined1")
  else {
    stop("At least one error term should be >0...")
  }
  if(is.null(type)) {
    type = allowedType[1]
  } else {
    if(! type %in% allowedType) stop("Your selected error model type ", type,
                                     " is not compatible with the given parameters. Allowed types: ", allowedType)
  }

  ## How do we calculate the standard deviation for these error models?
  myFunctions <- list(
    constant=function(ipred) {ipred*0 + add},
    proportional=function(ipred) {abs(prop*ipred)},
    combined1=function(ipred) {add + abs(prop*ipred)},
    combined2=function(ipred) {sqrt( add^2 + (prop*ipred)^2 )}
  )
  sigma <- myFunctions[[type]]
  res <- function(ipred, obs) {obs-ipred} ## residual error
  wres <- function(ipred, obs) { res(ipred, obs) / sigma(ipred) } ## Weighted residual error
  ll <- function(ipred, obs, log=TRUE) {
    sd <- sigma(ipred)
    i <- sd != 0 ## ignore the ones where sd = 0
    if(!any(i)) return(0)
    dnorm(ipred[i], obs[i], sd=sd[i], log=log)
  }
  inv_res <- function(ipred, res) { ipred + res }
  inv_wres <- function(ipred, wres) { inv_res(ipred, wres * sigma(ipred)) }
  return(structure(
    list(
      var=var,
      type=type,
      sigma=sigma,
      res=res,
      wres=wres,
      inv_res=inv_res,
      inv_wres=inv_wres,
      ll=ll,
      add=add,
      prop=prop
    ),
    class = c("error_model")
  ))
}

#' @export
print.error_model <- function(x, ...) {
  cat(x$var, " :", x$type, " ( ")
  if(x$add > 0) cat("add=", x$add, " ")
  if(x$prop > 0) cat("prop=", x$prop, " ")
  cat(")\n")
  invisible(x)
}
