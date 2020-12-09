library(benchmarkme)

getReference <- function() {
  file <- file.path( system.file("benchmark", package="tdmore"),
                     paste0(Sys.info()['nodename'], ".benchmark"))
  if(!file.exists(file)) {
    message("Running reference benchmark...")
    ref = benchmarkme::bm_prog_escoufier(runs = 1, verbose=FALSE)
    message("DONE. 1 jiffie = ", ref$elapsed, " s")
    saveRDS(ref, file)
  } else {
    ref = readRDS(file)
    message("Using cached benchmark result (1 jiffie = ", signif(ref$elapsed, 3), "s)...")
  }
  ref
}
ref <- getReference()

#' @noRd
#' @param lower The runtime should be higher than `lower*reference`
#' @param upper The runtime should be lower than `upper*reference`
expect_runtime <- function(expr, file, lower=0, upper=1.2) {


  res <- system.time(force(expr), gcFirst=TRUE)
  relTime <- res["elapsed"] / ref$elapsed #X jiffies

  if(isTRUE(as.logical(Sys.getenv("R_COVR")))) {
    rlang::warn("R_COVR is active, skipping runtime tests...")
    return()
  }

  if (!file.exists(file)) {
    cat(relTime, file = file)
    testthat::skip(message = "Created reference value")
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
