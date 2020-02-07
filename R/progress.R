## This is used to create dplyr-compatible Progress objects
MIN_WAIT_TIME <- 3

DplyrProgressFacade <- R6::R6Class("Progress",
   public = list(
     proxy=NULL,
     ## TODO: do not know if all of the below is really required...
     #n = NULL,
     #i = 0,

     initialize = function(n, min_time = 0, proxy=NULL, ...) {
       if(!is.null(proxy)) self$proxy <- proxy
       self$proxy$init(n)
       #self$n <- n
       #self$min_time <- min_time
     },
     tick = function() {
       self$proxy$step()
       self
     },
     stop = function() {
       self$proxy$term()
       self
     },
     print = function(...) { return(invisible(self)) }
   )
)

to_dplyr_progress <- function(x) {
  if(isTRUE(x) || (is.character(x) && x == "text"))
    return( dplyr::progress_estimated(n=1, min_time=MIN_WAIT_TIME) )
  if((is.character(x) && x == "none") || is.null(x))
    return( dplyr::progress_estimated(1, min_time=Inf) )

  if(is.list(x) && setequal(names(x), c("init", "step", "term")) ) {
    ## We need to make sure that calls to dplyr::progress_estimated
    ## are relayed to the init(N), step() and term() functions
    p <- DplyrProgressFacade$new(n=1, proxy=x)
    return(p)
  }
  if(R6::is.R6(x) && inherits(x, "Progress"))
    return(p)

  stop("Cannot translate .progress argument to a dplyr-compatible progress function...")
}
