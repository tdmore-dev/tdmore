#' Convert a numeric vector/list into a diagonal matrix.
#' If the vector/list is named, names will be used as the column and row names of the resulting diagonal matrix.
#'
#' @param vector a numeric vector/list
#'
#' @return a diagonal matrix
vectorToDiagonalMatrix <- function(vector) {
  if(length(vector)==1) {
    diagMatrix <- matrix(data=vector, nrow=1, ncol=1)
  } else {
    diagMatrix <- diag(x = vector)
  }
  if (!is.null(names(vector))) {
    colnames(diagMatrix) <- names(vector)
    rownames(diagMatrix) <- names(vector)
  }
  return(diagMatrix)
}

#' Melt prediction results.
#'
#' @param x dataframe returned by the predict function
#' @param se if FALSE, melts the data.frame as expected
#' if TRUE, assumes every output has 4 columns in the data.frame: XXX1, XXX1.median, XXX1.lower, XXX1.upper
#' And melts these into a variable (with values XXX1, XXX2, ...) and 4 columns: value, value.median, value.lower and value.upper
#'
#' @return the melted dataframe
#' @keywords internal
meltPredictions <- function(x, se=FALSE) {
  gathered <- x %>%
    tidyr::gather(key = "variable", value="value", -.data$TIME)
  if(!se) return(gathered)

  gathered %>% tidyr::extract(.data$variable, into=c("variable", "suffix"),
                              regex="^(.*?)(|\\.lower|\\.upper|\\.median)$") %>%
    dplyr::mutate(suffix=paste0("value", .data$suffix) ) %>%
    tidyr::spread(key=.data$suffix, value=.data$value)
}

#'
#' This function will check if the given package is loaded from a non-standard path.
#' If so, we are probably running in a sub-process of a test.
#'
#' If it is available, but in a standard path, then this is probably an outdated version and we need to install the current version.
#' If it is available, but as a dev_package, then we need to install the current version.
#'
#' The function always returns a temporary lib directory.
#'
#' @param pkgName the dev package name
#' @return a directory that can be added to libPaths. The directory is removed when the R session exits, or when tmp_lib is garbage-collected
#'
#' @export
ensurePackagePresent <- function(pkgName="tdmore") {
  tmp_lib <- tempfile("R_LIBS")
  dir.create(tmp_lib)
  reg.finalizer(tmp_lib, function(tmp_lib){
    unlink(tmp_lib, recursive = TRUE)
  }, onexit=TRUE)

  ## Is the library available on the current search path? If so, it probably was installed by covr or shinytest or ...
  ## If it is in R_LIBS_USER, then we refuse: this was installed through Build and Install, and is an old version!
  ## If no good version is available, then install the library in a temporary path and return this path
  pkgPat <- path.package(pkgName, quiet=T)
  if(is.null(pkgPat) || pkgName %in% devtools::dev_packages() || startsWith(pkgPat, Sys.getenv("R_LIBS_USER")) ) {
    #Either the package is not available
    #or it is currently in the dev packages
    #or we are about to load the (outdated) installed version!
    message("Installing temporary version of ", pkgName)
    pkg <- devtools::as.package(".")
    stopifnot(pkg$package == pkgName)
    utils::install.packages(repos = NULL,
                            lib = tmp_lib,
                            pkg$path,
                            type = "source",
                            INSTALL_opts = c("--example",
                                             "--install-tests",
                                             "--with-keep.source",
                                             "--with-keep.parse.data",
                                             "--no-multiarch"),
                            quiet = T)
  }
  return(tmp_lib)
}
