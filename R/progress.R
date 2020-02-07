## This is used to create dplyr-compatible Progress objects
MIN_WAIT_TIME <- 3

to_dplyr_progress <- function(x) {
  if(isTRUE(x) || (is.character(x) && x == "text"))
    return( dplyr::progress_estimated(n=1, min_time=MIN_WAIT_TIME) )
  if((isFALSE(x) || (is.character(x) && x == "none")) || (is.null(x)))
    return( dplyr::progress_estimated(1, min_time=Inf) )
  if(R6::is.R6(x) && inherits(x, "Progress"))
    return(x)

  stop("Cannot translate .progress argument to a dplyr-compatible progress function...")
}
