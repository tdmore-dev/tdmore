library(benchmarkme)

message("Running reference benchmark...")
ref = benchmarkme::bm_prog_escoufier(runs = 1, verbose=FALSE)
message("DONE. Reference elapsed time: ", ref$elapsed, " s")

#'
#' @param lower The runtime should be higher than `lower*reference`
#' @param upper The runtime should be lower than `upper*reference`
expect_runtime <- function(expr, file, lower=0, upper=1.2) {
  res <- system.time(force(expr), gcFirst=TRUE)
  relTime <- res["elapsed"] / ref$elapsed #X jiffies

  if (!file.exists(file)) {
    warning("Creating reference value", call. = FALSE)
    cat(relTime, file = file)
    testthat::succeed()
  }
  else {
    ref_val <- as.numeric( readChar(file, 99999) )
    testthat::expect_lt(relTime, ref_val*upper,
                               label=sprintf("Runtime for %s is longer than expected... %.2f jiffies (as compared to %.2f expected jiffies)", file, relTime, ref_val)
                                 )
    testthat::expect_gt(relTime, ref_val*lower,
                        label=sprintf("Runtime for %s is faster than expected... %.2f jiffies (as compared to %.2f expected jiffies)", file, relTime, ref_val)
    )
  }
}
