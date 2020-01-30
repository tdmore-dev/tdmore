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
#' @param ... extra arguments passed to stats::uniroot
#'
#' @return a recommendation object
#' @export
#' @importFrom stats uniroot
findDose <- function(tdmorefit, regimen=tdmorefit$regimen, doseRows=NULL, interval=c(0, 1E10), target, se.fit = FALSE, level = 0.95, mc.maxpts = 100, ...) {
  # Check if IOV is present in model
  iov <- tdmorefit$tdmore$iov
  if (!is.null(iov)) {
    #tdmorefitRegimen <- tdmorefit$regimen
    ## this is not required!
    #assert_that(getMaxOccasion(tdmorefitRegimen)  getMaxOccasion(regimen),
    #            msg="Number of occasions is different in tdmorefit regimen and findDose regimen")
  }
  if(is.null(doseRows))
    doseRows <- nrow(regimen)

  if(nrow(target) > 1) {
    stop("Cannot find the dose to hit multiple targets. Split the treatment regimen and perform findDose separately per target.")
  }
  if(ncol(target[, colnames(target) != "TIME", drop=FALSE]) > 1) {
    stop("Cannot find the dose to hit multiple targets. Split the treatment regimen and perform findDose separately per target.")
  }

  if (!se.fit) {
    # Find the best dose for the estimated parameters
    rootFunction <- function(AMT) {
      myRegimen <- updateRegimen(regimen = regimen, doseRows = doseRows, newDose = AMT)
      obs <- stats::predict(tdmorefit, newdata = target, regimen = myRegimen)
      result <- obs[, colnames(obs) != "TIME", drop=TRUE] - target[, colnames(target) != "TIME", drop=TRUE]
      if(length(result) > 1) stop("Cannot use findDose to hit multiple targets!")
      result
    }
    result <- runUniroot(rootFunction, interval, ...)
    return(convertResultToRecommendation(tdmorefit, result, regimen, doseRows, target))

  } else {
    # Find the dose for each Monte-Carlo sample
    mc <- sampleMC(tdmorefit, mc.maxpts = mc.maxpts)
    uniqueColnames <- make.unique(colnames(mc)) # needed for dplyr to have unique colnames

    mcResult <- purrr::map_dfr(mc$sample, function(i) {
      row <- mc[i,,drop=F] #make vector
      res <- unlist(row[-1]) # Remove 'sample' column
      names(res) <- names(coef(tdmorefit))

      mcRootFunction <- function(AMT) {
        myRegimen <- updateRegimen(regimen = regimen, doseRows = doseRows, newDose = AMT)
        obs <- predict.tdmore(object = tdmorefit$tdmore, newdata = target, regimen = myRegimen, parameters = res, covariates = tdmorefit$covariates)
        obs[, colnames(obs) != "TIME", drop=TRUE] - target[, colnames(target) != "TIME", drop=TRUE]
      }

      result <- runUniroot(mcRootFunction, interval, ...)
      names(row) <- uniqueColnames
      cbind(row, dose = result$root, f.root = result$f.root, iter = result$iter, estim.prec = result$estim.prec)
    })

    colnames(mcResult)[seq_len(length(uniqueColnames))] <- colnames(mc)

    return(convertMCResultToRecommendation(tdmorefit, mcResult, regimen, doseRows, target, level))
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
convertResultToRecommendation <- function(tdmorefit, result, regimen, doseRows, target) {
  return(structure(
    list(
      tdmorefit=tdmorefit,
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
convertMCResultToRecommendation <- function(tdmorefit, mcResult, regimen, doseRows, target, level) {
  ciLevel <- (1-level)/2
  dose <- c(
    dose.median = as.numeric(median(mcResult$dose)),
    dose.lower = as.numeric(quantile(mcResult$dose, ciLevel)),
    dose.upper = as.numeric(quantile(mcResult$dose, 1 - ciLevel))
  )
  return(structure(
    list(
      tdmorefit=tdmorefit,
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


#' Automatically guesses the target troughs for a given regimen
#'
#' @param regimen a treatment regimen
#' @param deltamin how much can we move the trough back to match an existing treatment time? in percentage of interdose-interval
#' @param deltaplus how much can we move the trough forward to match an existing treatment time? in percentage of interdose-interval
#' This method emits an error if any other treatment is found in the interval between `time` and `time+ii*(1-deltamin)`
#'
#' @return a numeric vector with the corresponding troughs
#' @export
getTroughs <- function(model, regimen, deltamin=1/4, deltaplus=1/4) {
  stopifnot( "FORM" %in% colnames(regimen) )
  getII <- function(x) {
    form <- tdmore::getMetadataByName(model, x)
    if(is.null(form)) stop("Formulation `", x, "' is not defined in the model metadata")
    form$dosing_interval
  }
  regimen$II <- purrr::map_dbl(regimen$FORM, getII)

  trough <- purrr::pmap_dbl(list(
    regimen$TIME,
    regimen$TIME + regimen$II*(1-deltamin),
    regimen$TIME + regimen$II,
    regimen$TIME + regimen$II*(1+deltaplus)),
    function(start, min, val, max) {
      x <- regimen$TIME
      error <- x > start & x < min #any treatments in the no-go zone?
      if(any(error)) stop("A treatment was detected between ", start, " and ", min, " at ", regimen$TIME[error])

      i <- x >= min & x <= max #is any existing treatment within the min-max interval?
      if(any(i)) {
        #if yes, use that treatment time as trough
        x[which.max(i)]
      } else {
        #if not, use the calculated trough
        val
      }
    }
  )
  trough
}

#' Optimize the regimen, using the metadata defined by the model
#'
#' This performs a step-wise optimization.
#'
#' @param fit tdmorefit object
#' @param regimen the treatment regimen to optimize
#'
#' @export
optimize <- function(fit, regimen=fit$regimen) {
  target <- list(
    TIME=getTroughs(fit$tdmore, regimen[regimen$FIX==FALSE, ])
  )
  targetMetadata <- tdmore::getMetadataByClass(fit$tdmore, "tdmore_target")
  outputVar <- fit$tdmore$res_var[[1]]$var
  target[outputVar] <- mean( c(targetMetadata$min, targetMetadata$max) )
  target <- tibble::as_tibble(target)

  #step-wise: row per row
  for(i in seq_along(target$TIME)) {
    row <- target[i, ]
    iterationRows <- which( regimen$FIX == FALSE & regimen$TIME < row$TIME )
    if(length(iterationRows) == 0) next
    rec <- findDose(fit, regimen, iterationRows, target=row)
    regimen <- rec$regimen
    regimen$FIX[iterationRows] <- TRUE #fix these rows in place!
  }

  return(structure(
    list(
      tdmorefit=tdmorefit,
      regimen = regimen,
      target=target
    ),
    class = c("recommendation")
  ))
}
