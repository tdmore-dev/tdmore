#' These functions perform estimate for multiple individuals at the same time.
#'
#' @param x a tbl with one or multiple rows
#' Each column from the tbl will be mapped to the `estimate()` function.
#' You can generate such a tbl using \link{data} functions.
#' @param ... additional arguments to add to each `estimate` call. If an argument
#' is present both here and in x, the argument here takes precedence.
#' @param .fit Either a string or NULL. If a string, the output will contain a column
#' with that name, storing the fit
#' @param .prediction Either a string or NULL. If a string, the output will contain a column
#' with that name, storing the predicted values for that fit
#' @param .elapsed Either a string or NULL. If a string, the output will contain a column
#' with that name, storing the required computation time (`System.time`) for that fit
#' @param .id column with the subject ID (to be ignored when calling `estimate`)
#'
#' @return a tibble with the requested columns for each fit
#'
#' @export
#' @examples
#' require(magrittr)
#' require(ggplot2)
#' data <- dataTibble(
#'     object=meropenem_model %>% tdmore(),
#'     observed=data.frame(ID=c(1,1,2,2), TIME=c(8,10,8,10), CONC=c(1.5,0.3,0.8,0.2))
#' )
#' fit <- data %>%
#'   posthoc(regimen=data.frame(TIME=0, AMT=1000))
#' fit %>% tidyr::unnest(c(observed, ipred), names_sep=".") %>% ggplot(aes(x=ipred.TIME)) +
#'   geom_point(aes(y=observed.CONC)) +
#'   geom_line(aes(y=ipred.CONC, group=ID)) +
#'   facet_wrap(~ID)
posthoc <- function(x, ..., .fit="fit", .prediction="ipred", .elapsed="elapsed", .id="ID") {
  if(!setequal(class(x), "list")) {
    res <- x %>% dplyr::rowwise() %>%
      dplyr::do({posthoc(., ..., .fit=.fit, .prediction=.prediction, .elapsed=.elapsed, .id=.id)})
    return(res)
  } else {
    stopifnot(is.list( x ))
    args <-x
    args[.id] <- NULL
    args[names(list(...))] <- list(...) #add/override using common arguments
    tic <- proc.time()
    fit <- do.call(estimate, c(args))
    toc <- proc.time()
    #Ensure a data.frame in the original input is spread over multiple columns...
    #You can try this yourself:
    #> foo <- list(bla1 = 1, bla2 = data.frame(TIME=3, CONC=2))
    #> as_tibble(foo)
    res <- lapply(x, function(x){if(rlang::is_atomic(x)) x else list(x)})
    if(!is.null(.fit)) res[[.fit]] <- list( fit )
    if(!is.null(.elapsed)) res[[.elapsed]] <- list( toc - tic )
    if(!is.null(.prediction)) res[[.prediction]] <- list( stats::predict(fit) )
    return( tibble::as_tibble(res) )
  }
}

#' @rdname posthoc
#' @param .obs Either a string or NULL. If a string, the output will contain a column
#' with that name, describing how many observations were taken into account
#' @section Proseval:
#' The proseval tool calculates a prospective evaluation, and therefore creates multiple outputs per input row.
#' If N is the number of observations, it will calculate N+1 fits.
#' The first fit will be the population fit (not using any observations).
#' The next fit will only use one observation.
#' The next fit will use two observations.
#' And so on.
#'
#' @export
proseval <- function(x, ..., .fit="fit", .prediction="ipred", .elapsed="elapsed", .obs="OBS") {
  if(is.null(.fit)) stop(".fit column needs to be present for proseval to work")
  if(tibble::is_tibble(x)) {
    res <- dplyr::rowwise(x) %>%
      dplyr::do({proseval(., ..., .fit=.fit, .prediction=.prediction, .elapsed=.elapsed, .obs=.obs)})
    return(res)
  } else {
    stopifnot(is.list( x ))

    output <- list()
    for(i in c(0, seq_len(nrow(x$observed))) ){
      argsPosthoc <- x
      if(i == 0) argsPosthoc['observed'] <- list(NULL)
      else argsPosthoc$observed <- x$observed[ seq_len(i), ]
      res <- posthoc(argsPosthoc, ..., .fit=.fit, .prediction=NULL, .elapsed=.elapsed)
      fit <- res[[.fit]][[1]]
      res[["observed"]] <- list(x$observed) #always use the same OBSERVED data.frame
      if(!is.null(.prediction)) res[[.prediction]] <- list( stats::predict(fit, x$observed) )
      if(!is.null(.obs)) res[[.obs]] <- i
      output[[i+1]] <- res
    }
    return(bind_rows(output))
  }
}

#' @rdname posthoc
#' @param optimize an optimization function for dose simulation, corresponding to function(fit, regimen, truth), and returning a list(nextTime, regimen)
#' @param predict a prediction function for dose simulation, corresponding to function(truth, newRegimen, nextTime)
#' If missing, we simply predict using the true fit.
#' We evaluate the nextObservations time points using the 'true' model and the new regimen.
#' @param .iterationFit column name for the iteration fit
#' @param .next_observed What will be the predicted next observation?
#' @param .next_regimen What will be the next regimen?
#' @param .verbose Be verbose when simulating?
#'
#' @section DoseSimulation:
#' This tool simulates a dose adaptation run for an individual.
#' We expect a .fit column in the data, that represents the best model fit for the given subject.
#' All other columns are passed to the posthoc function when calculating the next fit, except for
#' the 'observed' data (that is replaced by the simulated data)
#' and the 'regimen' data (that is replaced by the adapted regimen)
#'
#' The DoseSimulation iterates as follows:
#' 1) A fit is generated using the up-to-now observed concentrations
#' 2) The optimize method is called with the following arguments:
#'   optimize( fit, regimen, truth )
#'   It is up to the optimizer to determine how to optimize the treatment regimen.
#'   The optimizer returns a list( nextTime=nextObservations, regimen=newRegimen, extra=tibble() )
#'   If nextTime is not a finite number, we stop the simulation.
#'   If nextTime is missing, we use the next observation from the posthocfit.
#' 3) The predict method is called with the following arguments:
#'   predict(truth, newRegimen, nextTime)
#' This predicts the next observation.
#'
#'
#' @export
doseSimulation <- function(x, ..., optimize, predict,
                           .fit="fit", .iterationFit="iterationFit", .next_observed="next_observed", .next_regimen="next_regimen",
                           .verbose=TRUE) {
  if(! .fit %in% names(x)) stop(".fit column needs to be present in x for doseSimulation to work")
  if(missing(predict)) predict <- function(truth, regimen, time) {
    result <- stats::predict(truth, regimen=regimen, newdata=time)
    residuals.tdmorefit(truth, data=result, weighted=TRUE)
  }
  if(tibble::is_tibble(x)) {
    res <- dplyr::rowwise(x) %>%
      dplyr::do({doseSimulation(., ..., optimize=optimize, predict=predict,
                                .fit=.fit, .iterationFit=.iterationFit, .next_observed=.next_observed, .next_regimen=.next_regimen)})
    return(res)
  } else {
    stopifnot(is.list( x ))
    output <- list()

    truth <- x[[.fit]]
    regimen <- truth$regimen
    observed <- tibble::tibble(TIME=numeric()) ## specifying NULL will reuse the tdmorefit observed values!!
    repeat {
      if(.verbose) cat("Iteration #", nrow(observed), "\n")
      # 1) a fit is generated using up-to-now observed concentrations
      args <- x
      args[[.fit]] <- NULL
      args$observed <- observed
      args$regimen <- regimen
      res <- posthoc(args, ..., .fit=.iterationFit, .prediction=NULL, .elapsed=NULL)

      optimizeResult <- optimize(res[[.iterationFit]][[1]], regimen, truth)
      if(is.null(optimizeResult$nextTime)) optimizeResult$nextTime <- truth$observed$TIME[ nrow(observed)+1 ]

      extra <- optimizeResult$extra
      if(!is.null(extra)) res[ , names(extra) ] <- extra

      if(!all( is.finite(optimizeResult$nextTime) )) {
        output <- c(output, list(res) )
        break
      } else {
        regimen <- optimizeResult$regimen
        nextObservation <- predict(truth, regimen=regimen, time=optimizeResult$nextTime)
        observed <- bind_rows(observed,
                              model.frame.tdmorefit(truth, data=nextObservation)
        )
        res[[.next_observed]] <- list(observed)
        res[[.next_regimen]] <- list(regimen)
        output <- c(output, list(res) )
      }
    }

    outTibble <- bind_rows(output)

    res <- lapply(x, function(x){if(rlang::is_atomic(x)) x else list(x)})
    outTibble[, names(res)] <- res

    outTibble
  }
}
