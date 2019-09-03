## Optimize the treatment regimen in funcion of the given target

#' Plot the fitted prediction under all discrete regimen possibilities
#'
#' @param tdmorefit A tdmore fit
#' @param regimen Candidate regimen
#' @param searchspace a list with searchSpaces for every row in the regimen table, or a single searchspace
#' In case this is a single searchSpace, it will be attributed to
#' every row in the regimen table after the last observation.
#' Any dose rows where ADDL is spanning over the last observation will
#' not be attributed to the searchSpace.
#' @param evaluators named list of evaluation functions(tdmorefit, regimen)
#'
#' @importFrom ggplot2 ggplot geom_line geom_hline aes_
#' @importFrom gridExtra tableGrob
#'
#' @return list with source data for plots, as well as the plots in question
#' @export
#'
plotpossibilities <- function(tdmorefit, regimen, searchspace, evaluators) {
  x <- evaluate(tdmorefit, regimen, searchspace, evaluators)

  db <- x$pred
  eval <-x$eval

  p1 <- ggplot(db) +
    geom_line(aes_(x="TIME", y="CONC")) +
    geom_hline(yintercept=0.05)

  p2 <- tableGrob(eval)

  list(pkPlot=p1, tablePlot=p2, pred=db, table=eval)
}

#' Evaluate all discrete regimen possibilities
#'
#' @param tdmorefit A tdmore fit
#' @param regimen Candidate regimen
#' @param searchspace a list with searchSpaces for every row in the regimen table, or a single searchspace
#' In case this is a single searchSpace, it will be attributed to
#' every row in the regimen table after the last observation.
#' Any dose rows where ADDL is spanning over the last observation will
#' not be attributed to the searchSpace.
#' @param evaluators named list of evaluation functions(tdmorefit, regimen)
#'
#' @return list with evaluation metrics per scenario, and model predictions
#' @export
#'
evaluate <- function(tdmorefit, regimen, searchspace, evaluators) {
  ## In case we only provide a single searchSpace, attribute it to all treatment rows occuring after
  if(inherits(searchspace, "searchspace")) {
    tmax <- max(tdmorefit$observed$TIME)

    ## Split any SS regimens
    lastDose <- regimen$TIME + regimen$II*regimen$ADDL

    tmp <- list()
    # for(i in seq_len(nrow(regimen))) {
    #   tmp[[i]] <-
    # }
    tmp[[ which(lastDose > tmax) ]] <- searchspace
    searchspace <- tmp
  }
  sapply(searchspace, function(x) {assert_that(is.null(x) || inherits(x, "searchspace"))})

  ## We need a table with all possibilities for all rows
  ## 1.AMT 1.II 2.AMT 3.AMT
  ## O1   O1    O1    O1
  ## O1   O1    O1    O2
  ## O1   O1    O2    O1
  ## etc.

  grid <- possibilitygrid(searchspace)

  ## Main loop
  newGrid <- apply(grid, 1, function(x) {
    myRegimen <- transformRegimen(regimen, x)
    for(j in names(evaluators)) {
      evaluatorFunction <- evaluators[[j]]
      x[j] <- evaluatorFunction(tdmorefit, myRegimen)
    }
    print(myRegimen)
    list(eval=x, pred=predict.tdmorefit(tdmorefit, regimen=myRegimen))
  })
  names(newGrid) <- seq_along(newGrid)

  db <- purrr::map_dfr(newGrid, function(x) {x$pred})
  eval <- purrr::map_dfr(newGrid, function(x) {x$eval})

  list(pred=db, table=eval)
}

#' For a given set of planned times, round to the nearest reference time.
#'
#' This can be used to move an intended administration to fit the historic treatment regimen of a
#' patient in retrospective data.
#'
#' @param x a vector of times to round
#' @param reference a vector of times, used to round
#' @param range how much can be subtracted or added from/to 'x' to get to the reference
#' @param if.na if no corresponding reference time is found, what should we do?
#' if 'keep', then the original time is kept
#' if 'NA', then we set the value to NA
#'
#' @return A vector of the same length as x. Values close to reference values are replaced by the reference value.
#'
#' @export
#'
roundTime <- function(x, reference, range=c(-0.5, 0.5), if.na=c("keep", "NA")) {
  if.na <- match.arg(if.na)
  lower <- range[1]
  upper <- range[2]
  for(i in seq_along(x)) {
    delta <- reference - x[i]
    ## if reference is higher, then delta is positive
    ## so delta should be smaller than 'upper'

    matches <- delta >= lower & delta <= upper
    if(!any(matches)) {
      # no match
      if(if.na=="NA") x[i] <- NA
    } else {
      # match !
      delta[!matches] <- Inf
      j <- which.min( abs(delta) ) #which one is closest?
      x[i] <- reference[j]
    }
  }
  x
}
