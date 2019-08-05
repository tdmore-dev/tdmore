#' This function performs a posthoc estimate for multiple individuals at the same time.
#' All parameters may have an additional ID column to identify every individual.
#'
#' @inheritParams estimate
#' @param ... passed to \link{estimate} function
#' @param .prediction Either a string or NULL. If a string, the output will contain a column
#' with that name, storing the predicted values for that specific fit.
#'
#' @return a tibble with columns ID,
#' and column `fit` with the tdmorefit object
#'
#' @export
posthoc <- function(object, observed=NULL, regimen=NULL, covariates=NULL, par=NULL, fix=NULL, ..., .prediction="ipred") {
  ._posthoc(object, observed=observed, regimen=regimen, covariates=covariates, par=par, fix=fix, ..., .proseval=FALSE, .prediction=.prediction)
}

#' This function calculates a prospective evaluation.
#' For every ID, it will calculate as many fits as there are rows of observations +1.
#' The first fit will be the population fit (not using any observations).
#' The next fit will only use one observation.
#' The next fit will use an extra observation.
#' And so on.
#'
#' @inheritParams posthoc
#'
#' @return a tibble with columns ID,
#' OBS (showing how many observations were included in this fit),
#' and column `fit` with the tdmorefit object (nested in a list)
#'
#' @export
proseval <- function(object, observed=NULL, regimen=NULL, covariates=NULL, par=NULL, fix=NULL, ..., .prediction="ipred") {
  ._posthoc(object, observed=observed, regimen=regimen, covariates=covariates, par=par, fix=fix, ..., .proseval=TRUE, .prediction=.prediction)
}

._posthoc <- function(object, observed=NULL, regimen=NULL, covariates=NULL, par=NULL, fix=NULL, ..., .proseval=FALSE, .prediction=NULL) {
  check <- function(x) {
    assert_that(
      is.null(x) || {
        is.data.frame(x)
      }
    )
  }
  arguments <- list(observed=observed, regimen=regimen, covariates=covariates, par=par, fix=fix)
  purrr::map(arguments, check)

  getId <- function(.x) {
    if(is.null(.x)) {c()} else {.x$ID}
  }
  ID <- arguments %>% purrr::map(getId) %>% unlist %>% unique

  getForId <- function(.x, id) {
    if(is.null(.x)) return(NULL)
    if(! "ID" %in% colnames(.x)) return(.x)
    .x %>% dplyr::filter(.data$ID == id) %>% dplyr::select(-.data$ID)
  }
  doEstimate <- function(argumentsPerId) {
    fit <- do.call(estimate, args=c(list(object), argumentsPerId, list(...)))
    do.call( tibble::tibble, args=c(
      #actually, all of these arguments are already in the tdmorefit object!
      #map( argumentsPerId, ~list(.x)),
      list(fit=list(fit), .rows=1)
    ))
  }

  pb <- dplyr::progress_estimated(n=length(ID))
  estimateForId <- function(data, key){
    pb$tick()$print()
    argumentsPerId <- arguments %>% purrr::map(~ getForId(.x, key$ID))
    if(.proseval) {
      if(is.null(argumentsPerId$observed) ||
         !is.data.frame(argumentsPerId$observed) ||
         !"TIME" %in% colnames(argumentsPerId$observed))
        stop("Please specify `observed` to use proseval")

      result <- list()
      res <- NULL
      for(i in c(0, seq_along(argumentsPerId$observed$TIME)) ) {
        iterationArguments <- argumentsPerId
        iterationArguments$observed <- if(i==0) NULL else iterationArguments$observed[1:i, ]
        ## Start estimation at result from last time
        if(!is.null(res) && is.null(iterationArguments$par)) {
          iterationArguments$par <- coef(res$fit[[1]])
        }
        res <- doEstimate(iterationArguments)
        res$OBS <- i
        if(!is.null(.prediction))
          res[,.prediction] <- list(list( predict(res$fit[[1]], argumentsPerId$observed) ))

        result[[i+1]] <- res
      }
      dplyr::bind_rows(result)
    } else {
      res <- doEstimate(argumentsPerId)
      if(!is.null(.prediction))
        res[,.prediction] <- list(list( predict(res$fit[[1]], argumentsPerId$observed) ))
      res
    }
  }
  res <- tibble::tibble(ID=ID) %>% dplyr::group_by(ID) %>% dplyr::group_modify(estimateForId)
  res
}
