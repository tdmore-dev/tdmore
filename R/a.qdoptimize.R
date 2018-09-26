#Quick and Dirty optimize
#' Find the dose required
#'
#' @param tdmorefit the tdmorefit object
#' @param regimen The regimen
#' @param doseRows Which rows of the regimen to adapt when searching for a new dose, or NULL to take the last one
#' @param interval Which interval to search a dose in. Defaults to a ridiculously high range
#' @param target Target value, as a data.frame
#' @param ... Extra arguments passed to stats::uniroot
#'
#' @return A list with at least four components:
#' \code{root} and \code{f.root} give the location of the root and the value of the function evaluated
#' at that point. \code{iter} and \code{estim.prec} give the number of iterations used
#' and an approximate estimated precision for \code{root}.
#'  (If the root occurs at one of the endpoints, the estimated precision is \code{NA}.)
#' @export
#' @importFrom stats uniroot
findDose <- function(tdmorefit, regimen, doseRows=NULL, interval=c(0, 1E10), target, ...) {
  if(is.null(doseRows)) doseRows <- nrow(regimen)

  dose <- numeric(length = length(doseRows))
  names(dose) <- paste0(doseRows, ".AMT")

  rootFunction <- function(AMT) {
    dose <- dose*0 + AMT
    myRegimen <- transform(regimen, dose)
    obs <- predict(tdmorefit, newdata=target, regimen=myRegimen)
    obs[ , colnames(obs) != "TIME"] - target[, colnames(target) != "TIME"]
  }
  iValues <- c(rootFunction(interval[1]), rootFunction(interval[2]) )
  if(sign(iValues[1]) == sign(iValues[2])) {
    warning("Predicted values at edges of interval both ",
            switch(sign(iValues[1]), "0"="equal to", "1"="above", "2"="below"),
            " target, returning closest value...")
    i <- which.min(abs(iValues))
    return(list(
      root=interval[i],
      f.root=iValues[i],
      iter=0,
      estim.prec=Inf
    ))
  }

  uniroot(f=rootFunction, interval=interval)
}
