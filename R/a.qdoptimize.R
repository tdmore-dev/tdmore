#' Find the dose required (quick and dirty, to be optimised)
#'
#' @param tdmorefit the tdmorefit object
#' @param regimen the regimen
#' @param doseRows which rows of the regimen to adapt when searching for a new dose, or NULL to take the last one
#' @param interval which interval to search a dose in. Defaults to a ridiculously high range
#' @param target target value, as a data.frame
#' @param se.fit TRUE to provide a confidence interval on the dose prediction, adding columns dose.median, dose.lower and dose.upper
#' @param level the confidence interval on the dose, only used if se.fit is true
#' @param mc.maxpts maximum number of points to sample in Monte Carlo simulation
#' @param .progress see plyr::ddply
#' @param .parallel see plyr::ddply
#' @param ... extra arguments passed to stats::uniroot
#'
#' @return a recommendation object
#' @export
#' @importFrom stats uniroot
findDose <- function(tdmorefit, regimen, doseRows=NULL, interval=c(0, 1E10), target, se.fit = FALSE, level = 0.95, mc.maxpts = 100, .progress="none", .parallel=FALSE, ...) {
    if(!se.fit) {
      # Find the best dose for the estimated parameters
      rootFunction <- function(AMT) {
        myRegimen <- updateRegimen(regimen = regimen, doseRows = doseRows, newDose = AMT)
        obs <- predict(tdmorefit, newdata=target, regimen=myRegimen)
        obs[ , colnames(obs) != "TIME"] - target[, colnames(target) != "TIME"]
      }
      result <- runUniroot(rootFunction, interval, ...)
      return(convertResultToRecommendation(result, regimen, doseRows, target))

    } else {
      # Find the dose for each Monte-Carlo sample
      mc <- as.data.frame( mnormt::rmnorm(mc.maxpts, mean=coef(tdmorefit), varcov=vcov(tdmorefit)) )
      mc$sample <- 1:mc.maxpts
      parameters <- coef(tdmorefit)

      mcResult <- plyr::ddply(mc, 1, function(row) {
        paramValues <- row[names(parameters)]

        mcRootFunction <- function(AMT) {
          myRegimen <- updateRegimen(regimen = regimen, doseRows = doseRows, newDose = AMT)
          obs <- predict.tdmore(object=tdmorefit$tdmore, newdata=target, regimen=myRegimen, parameters=unlist(paramValues), covariates=tdmorefit$covariates)
          obs[, colnames(obs) != "TIME"] - target[, colnames(target) != "TIME"]
        }

        result <- runUniroot(mcRootFunction, interval, ...)
        cbind(row, dose=result$root, f.root=result$f.root, iter=result$iter, estim.prec=result$estim.prec)
      }, .progress=.progress, .parallel=.parallel)

      return(convertMCResultToRecommendation(mcResult, regimen, doseRows, target, level))
    }
}

runUniroot <- function(rootFunction, interval, ...) {
  iValues <- c(rootFunction(interval[1]), rootFunction(interval[2]))

  if(sign(iValues[1]) == sign(iValues[2])) {
    warning("Predicted values at edges of interval both ",
            switch(sign(iValues[1]), "0"="equal to", "1"="above", "2"="below"),
            " target, returning closest value...")
    i <- which.min(abs(iValues))

    result <- list(
      root=interval[i],
      f.root=iValues[i],
      iter=0,
      estim.prec=Inf
    )
    return(result)
  } else {
    return(uniroot(f=rootFunction, interval=interval, ...))
  }
}

#' Convert the result of the uniroot function into a recommendation object.
#'
#' @param result the uniroot result
#' @param regimen the regimen
#' @param doseRows which rows of the regimen to adapt when searching for a new dose, or NULL to take the last one
#' @param target the original target
#'
#' @return a recommendation object
#' @keywords internal
convertResultToRecommendation <- function(result, regimen, doseRows, target) {
  return(structure(
    list(
      dose = result$root,
      regimen = updateRegimen(
        regimen = regimen,
        doseRows = doseRows,
        newDose = result$root
      ),
      target=target,
      result = result
    ),
    class = c("recommendation")
  ))
}

#' Convert the Monte-Carlo uniroot results into a recommendation object.
#'
#' @param mcResult a dataframe with all the uniroot results, one per Monte-Carlo sample
#' @param regimen the regimen
#' @param doseRows which rows of the regimen to adapt when searching for a new dose, or NULL to take the last one
#' @param target the original target
#' @param level the confidence interval on the dose
#'
#' @return a recommendation object
#' @importFrom dplyr summarise
#' @keywords internal
convertMCResultToRecommendation <- function(mcResult, regimen, doseRows, target, level) {
  ciLevel <- (1-level)/2
  dose <- c(
    dose.median = as.numeric(median(mcResult$dose)),
    dose.lower = as.numeric(quantile(mcResult$dose, ciLevel)),
    dose.upper = as.numeric(quantile(mcResult$dose, 1 - ciLevel))
  )
  return(structure(
    list(
      dose = dose,
      regimen = updateRegimen(
        regimen = regimen,
        doseRows = doseRows,
        newDose = dose['dose.median']
      ),
      result = mcResult,
      target=target
    ),
    class = c("recommendation")
  ))
}

#' @export
as.double.recommendation <- function(x, ...) {
  x$dose
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

  dose <- numeric(length = length(doseRows)) + as.numeric(newDose)
  names(dose) <- paste0(doseRows, ".AMT")
  updatedRegimen <- transformRegimen(regimen, dose)

  return(updatedRegimen)
}
