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
  if(class(searchspace) == "searchspace") {
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
  sapply(searchspace, function(x) {assert_that(is.null(x) || class(x) == "searchspace")})

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
