#' A facade object to use Shiny Progress as if they are a progress::progress_bar
#'
#' @description
#' Used to translate dplyr progress to a Shiny progress bar
#' @noRd
ShinyToProgressFacade <- R6::R6Class("progress_bar",
  public = list(
    #' @field proxy Shiny Progress object
    proxy=NULL,

    #' @field total total ticks for progress; can be reinitialized
    total = NULL,

    #' @field i current tick for progress; can be reinitialized
    i = 0,

    #' @field offset offset for Shiny Progress
    offset=0,

    #' @description
    #' Initialize a shiny progress to dplyr facade. This can either be called with the `proxy` and `amount`
    #' arguments for a first initialization (e.g. when using `new`), or with
    #' `n` and `min_time` arguments for the dplyr-like initialization
    #'
    #' @param total Total number of items
    #' @param min_time Progress bar will wait until at least `min_time`
    #'   seconds have elapsed before displaying any results.
    #' @param proxy Shiny Progress object
    #' @param ... ignored
    #' @return an object that can be used in `.progress` arguments for dplyr
    initialize = function(total=1, proxy=NULL, ...) {
      if(!is.null(proxy)) {
        self$proxy <- proxy
        self$offset <- proxy$getValue()
        if(is.null(self$offset)) self$offset <- 0 #no value defined yet
      }
      self$total <- total
      self$i <- 0
    },

    #' @description
    #' Advance a tick
    tick = function() {
      self$i = self$i + 1
      self$proxy$set(value=self$offset + self$i / self$total )
      invisible(self)
      self
    }
  )
)


## This is used to create dplyr-compatible Progress objects
MIN_WAIT_TIME <- 1

# R6 does not have any way to find out if an object is instance of a class
# This heuristic is as close as we can get: does the public description match?
matchesR6 <- function(x, class) {
  R6::is.R6(x) && inherits(x, class$classname) &&
    setequal(c( names(class$public_fields), names(class$public_methods), ".__enclos_env__" ), names(x)) #all public fields match
}

## This generates a `progress_bar` from the progress package
to_progress <- function(x) {
  if (!requireNamespace("progress", quietly = TRUE)) {
    warning("the `progress` package must be installed to show progress bars")
    p <- list(initialize=function(...){}, tick=function(...){})
    return(p)
  }

  if(isTRUE(x) || (is.character(x) && x == "text"))
    return(
      progress::progress_bar$new(total=1, show_after=MIN_WAIT_TIME)
      )
  if((isFALSE(x) || (is.character(x) && x == "none")) || (is.null(x)))
    #progress bar that is never shown
    # Inf is not supported by progress
    return( progress::progress_bar$new(total=1, show_after=1e12) )
  if(R6::is.R6(x) && inherits(x, "Progress")) {
    # it could be a dplyr::Progress object or a shiny::Progress object...
    if(matchesR6(x, shiny::Progress)) {
      return( ShinyToProgressFacade$new(proxy=x) )
    } else {
      stop("dplyr progress_estimated is deprecated")
    }
  }
  if(R6::is.R6(x) && inherits(x, "progress_bar")) {
    if(matchesR6(x, progress::progress_bar)) {
      return(x)
    } else {
      stop("unknown progress bar type progress_bar")
    }
  }

  stop("Cannot translate .progress argument to a progress-compatible progress function...")
}
