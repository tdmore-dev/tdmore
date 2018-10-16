#' Find the dose required (quick and dirty, to be optimised)
#'
#' @param tdmorefit the tdmorefit object
#' @param regimen the regimen
#' @param doseRows which rows of the regimen to adapt when searching for a new dose, or NULL to take the last one
#' @param interval which interval to search a dose in. Defaults to a ridiculously high range
#' @param target target value, as a data.frame
#' @param ... extra arguments passed to stats::uniroot
#'
#' @return a recommendation object
#' @export
#' @importFrom stats uniroot
findDose <- function(tdmorefit, regimen, doseRows=NULL, interval=c(0, 1E10), target, ...) {

  rootFunction <- function(AMT) {
    myRegimen <- updateRegimen(regimen = regimen, doseRows = doseRows, newDose = AMT)
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

  result <- uniroot(f=rootFunction, interval=interval)

  structure(list(
    dose = result$root,
    regimen = updateRegimen(
      regimen = regimen,
      doseRows = doseRows,
      newDose = result$root
    ),
    result = result
  ),

  class = c("recommendation"))
}

#' Update a regimen with the specified dose.
#'
#' @param regimen the regimen to update
#' @param doseRows which rows of the regimen to adapt, if not specified, the last dose will be adapted
#' @param newDose the specified new dose
#'
#' @return the updated regimen
#' @export
updateRegimen <- function(regimen, doseRows = NULL, newDose) {
  if (is.null(doseRows))
    doseRows <- nrow(regimen)

  dose <- numeric(length = length(doseRows)) + newDose
  names(dose) <- paste0(doseRows, ".AMT")
  updatedRegimen <- transform(regimen, dose)

  return(updatedRegimen)
}
