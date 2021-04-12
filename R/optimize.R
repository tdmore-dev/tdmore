#' @describeIn doseSimulation Find the dose required to hit a target.
#'
#' @param fit the tdmorefit object
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
findDose <- function(fit, regimen=fit$regimen, doseRows=NULL, interval=c(0, 1E10), target, se.fit = FALSE, level = 0.95, mc.maxpts = 100, ...) {
  # Check if IOV is present in model
  iov <- fit$tdmore$iov
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
      obs <- stats::predict(fit, newdata = target, regimen = myRegimen)
      result <- obs[, colnames(obs) != "TIME", drop=TRUE] - target[, colnames(target) != "TIME", drop=TRUE]
      if(length(result) > 1) stop("Cannot use findDose to hit multiple targets!")
      result
    }
    result <- runUniroot(rootFunction, interval, ...)
    return(convertResultToRecommendation(fit, result, regimen, doseRows, target))

  } else {
    # Find the dose for each Monte-Carlo sample
    mc <- sampleMC(fit, mc.maxpts = mc.maxpts)$mc

    result <- purrr::map(mc$sample, function(i) {
      row <- mc[i,,drop=F] #make vector
      res <- unlist(row[-1]) # Remove 'sample' column
      names(res) <- names(coef(fit))

      mcRootFunction <- function(AMT) {
        myRegimen <- updateRegimen(regimen = regimen, doseRows = doseRows, newDose = AMT)
        obs <- stats::predict(object = fit$tdmore, newdata = target, regimen = myRegimen, parameters = res, covariates = fit$covariates)
        obs[, colnames(obs) != "TIME", drop=TRUE] - target[, colnames(target) != "TIME", drop=TRUE]
      }

      runUniroot(mcRootFunction, interval, ...)
    })

    return(convertMCResultToRecommendation(fit, mc, result, regimen, doseRows, target, level))
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
#' @noRd
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
      result = list(result)
    ),
    class = c("recommendation")
  ))
}

#' Convert the Monte-Carlo uniroot results into a recommendation object.
#'
#' @param mc a dataframe with all the monte-carlo samples for parameters
#' @param result list with uniroot results
#' @param regimen the regimen
#' @param doseRows which rows of the regimen to adapt when searching for a new dose, or NULL to take the last one
#' @param target the original target
#' @param level the confidence interval on the dose
#'
#' @return a recommendation object
#' @importFrom dplyr summarise
#' @keywords internal
#' @noRd
convertMCResultToRecommendation <- function(tdmorefit, mc, result, regimen, doseRows, target, level) {
  doses <- purrr::map_dbl(result, ~.x$root)

  ciLevel <- (1-level)/2
  dose <- c(
    dose.median = as.numeric(median(doses)),
    dose.lower = as.numeric(quantile(doses, ciLevel)),
    dose.upper = as.numeric(quantile(doses, 1 - ciLevel))
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
      result=result,
      target=target
    ),
    class = c("recommendation")
  ))
}

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
#' @noRd
updateRegimen <- function(regimen, doseRows = NULL, newDose) {
  if (is.null(doseRows))
    doseRows <- nrow(regimen)

  dose <- numeric(length = length(doseRows)) + as.numeric(newDose)
  names(dose) <- paste0(doseRows, ".AMT")
  updatedRegimen <- transformRegimen(regimen, dose)

  return(updatedRegimen)
}

#' This changes only the mantissa in a double number
#' We can use this to ensure something is "before" another number,
#' while only having an infinitely small difference to that number.
#' @noRd
modifyMantissa <- function(x, a=-2^-52) {
  sign <- sign(x)
  x <- abs(x)
  exponent <- floor(log2(x))
  mantissa <- x * 2^(-exponent)

  x <- mantissa * 2^exponent

  mantissa2 <- mantissa + a #52 bits for mantissa
  x2 <- mantissa2 * 2^exponent
  x2[sign == 0] <- a
  x2[sign == -1] <- -1 * x2

  # x2 is 1 single bit lower than x
  # sprintf("%a", x)
  # sprintf("%a", x2)
  #

  x2
}

#' @describeIn doseSimulation Automatically guesses the target troughs for a given regimen
#'
#' @param model tdmore model with formulation metadata
#' @param regimen a treatment regimen
#' @param deltamin how much can we move the trough back to match an existing treatment time? in percentage of interdose-interval
#' @param deltaplus how much can we move the trough forward to match an existing treatment time? in percentage of interdose-interval
#' This method emits an error if any other treatment is found in the interval between `time` and `time+ii*(1-deltamin)`
#' @param adj some prediction models will display the peak rather than the trough if we use the exact treatment time.
#' To counteract this, we subtract an infitisemal small number. Specify an actual number,
#' or specify TRUE to subtract a single bit from the mantissa (the lowest theoretical amount to ensure the returned
#' value is lower than the treatment TIME)
#'
#' @return a numeric vector with the corresponding troughs
#' @engine
getTroughs <- function(model, regimen, deltamin=1/4, deltaplus=1/4, adj=TRUE) {
  stopifnot( "FORM" %in% colnames(regimen) )
  regimen$II <- getDosingInterval(regimen$FORM, model)

  trough <- purrr::pmap_dbl(list(
    regimen$TIME,
    regimen$TIME + regimen$II*(1-deltamin),
    regimen$TIME + regimen$II,
    regimen$TIME + regimen$II*(1+deltaplus)),
    function(start, min, val, max) {
      x <- regimen$TIME
      error <- x > start & x < min #any treatments in the no-go zone?
      if(any(error)) {
        warning("A treatment was planned/administered earlier than allowed by the inter-dose interval. ",
        "Treatment at ", regimen$TIME[error], " is planned in zone [", start, ", ", min, "]. We will use that treatment time as the 'trough'. ")
        return(x[ error ][1])
      }

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
  if(is.numeric(adj)) return( trough - adj )
  if(isTRUE(adj)) {
    # see https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html for background
    # https://stackoverflow.com/questions/50217954/double-precision-64-bit-representation-of-numeric-value-in-r-sign-exponent
    return( vapply(trough, modifyMantissa, FUN.VALUE=numeric(1), a=-2^-52) )
  }
  trough - adj
}

#' @describeIn doseSimulation Optimize the regimen, using the metadata defined by the model
#'
#' This performs a step-wise optimization of multiple doses.
#'
#' @param fit tdmorefit object
#' @param regimen the treatment regimen to optimize
#' @param targetMetadata defined target troughs as list(min=X, max=Y). If NULL or all NA, taken from the model metadata.
#'
#' @export
findDoses <- function(fit, regimen=fit$regimen, targetMetadata=NULL) {
  if(! "FIX" %in% colnames(regimen) ) regimen$FIX <- FALSE
  if(is.null(targetMetadata) || all(is.na(targetMetadata))) {
    targetMetadata <- tdmore::getMetadataByClass(fit$tdmore, "tdmore_target")
    if(is.null(targetMetadata)) stop("No target defined in model metadata")
  }
  stopifnot( all( c("min", "max") %in% names(targetMetadata) ) )
  #stopifnot( all( c("min", "max", "deltamin", "deltaplus") %in% names(targetMetadata) ) )

  # target <- list(
  #   TIME=getTroughs(fit$tdmore, regimen[regimen$FIX==FALSE, ], adj=TRUE, deltamin=targetMetadata$deltamin, deltaplus = targetMetadata$deltaplus)
  # )
  target <- list(
    TIME=getTroughs(fit$tdmore, regimen[regimen$FIX==FALSE, ], adj=TRUE)
  )


  outputVar <- fit$tdmore$res_var[[1]]$var
  targetValue <- mean( c(targetMetadata$min, targetMetadata$max) )
  if(is.na(targetValue)) stop("Target not defined, cannot optimize treatment...")
  target[outputVar] <- targetValue
  target <- tibble::as_tibble(target)

  modified <- rep(FALSE, nrow(regimen))
  result <- list()
  #step-wise: row per row
  for(i in seq_along(target$TIME)) {
    row <- target[i, ]
    iterationRows <- which( regimen$FIX == FALSE & regimen$TIME < row$TIME & !modified )
    if(length(iterationRows) == 0) next
    rec <- findDose(fit, regimen, iterationRows, target=row)
    regimen <- rec$regimen
    roundedAmt <- purrr::pmap_dbl(list(regimen$AMT, regimen$FORM), function(amt, form) {
      form <- tdmore::getMetadataByName(fit$tdmore, form)
      form$round_function(amt)
    })
    regimen$AMT[iterationRows] <- roundedAmt[iterationRows] #rounded amounts only
    modified[ iterationRows ] <- TRUE
    result <- c( result, rec$result )
  }

  return(structure(
    list(
      tdmorefit=fit,
      dose=regimen$AMT[ modified ],
      regimen = regimen,
      target=target,
      result=result
    ),
    class = c("recommendation")
  ))
}
