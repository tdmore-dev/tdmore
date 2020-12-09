#' Convert a numeric vector/list into a diagonal matrix.
#' If the vector/list is named, names will be used as the column and row names of the resulting diagonal matrix.
#'
#' @param vector a numeric vector/list
#'
#' @return a diagonal matrix
#' @noRd
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

#' This utility function will check if the given package is loaded from a non-standard path.
#' If so, we are probably running in a sub-process of a test.
#'
#' If it is available, but in a standard path, then this is probably an outdated version and we need to install the current version.
#' If it is available, but as a dev_package, then we need to install the current version.
#'
#' The function always returns a temporary lib directory.
#'
#' @param pkgName the dev package name
#' @param quiet FALSE to display log messages
#' @return a directory that can be added to libPaths
#'
#' @export
#' @keywords internal
ensurePackagePresent <- function(pkgName="tdmore", quiet=TRUE) {
  tmp_lib <- tempfile("R_LIBS")
  dir.create(tmp_lib)

  ## Is the library available on the current search path? If so, it probably was installed by covr or shinytest or ...
  ## If it is in R_LIBS_USER, then we refuse: this was installed through Build and Install, and is an old version!
  ## If no good version is available, then install the library in a temporary path and return this path
  pkgPat <- path.package(pkgName, quiet=T)
  r_libs_user <- Sys.getenv("R_LIBS_USER")
  if(nchar(r_libs_user) == 0) r_libs_user <- "DOES_NOT_EXIST"
  if(is.null(pkgPat) || pkgName %in% devtools::dev_packages() || startsWith(pkgPat, r_libs_user) ) {
    if(!quiet) message("Package ", pkgName, appendLF=F)
    if(!quiet && is.null(pkgPat)) message(" is not available. ", appendLF=F)
    else if(!quiet && pkgName %in% devtools::dev_packages()) message(" is loaded as dev package. ", appendLF=F)
    else if(!quiet && startsWith(pkgPat, r_libs_user)) message(" is available in R_LIBS_USER. ", appendLF=F)
    if(!quiet) message("Installing new version...")
    #Either the package is not available
    #or it is currently in the dev packages
    #or we are about to load the (outdated) installed version!
    message("Installing temporary version of ", pkgName)
    pkg <- devtools::as.package(".")
    stopifnot(pkg$package == pkgName)
    callr::r( ## utils::install.packages refuses if the package is already loaded
      func = function(...) {
        utils::install.packages(...) ## return function value does not work in callr
        return(NULL)
        },
       args=list(
         pkgs=pkg$path,
         repos = NULL,
              lib = tmp_lib,
              type = "source",
              INSTALL_opts = c(#"--example",
                               #"--install-tests",
                               "--with-keep.source",
                               "--with-keep.parse.data",
                               "--no-docs",
                               "--no-html",
                               "--no-lock",
                               "--no-help", "--no-demo", "--no-exec", "--data-compress=none",
                               "--no-byte-compile", "--no-staged-install",
                               "--no-test-load",
                               "--no-multiarch"),
              quiet = quiet)
    )
  } else {
    if(!quiet) message("Package ", pkgName, " is available in temporary library, not installing new version...")
  }
  return(tmp_lib)
}
