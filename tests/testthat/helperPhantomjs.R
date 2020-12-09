if(.Platform$OS.type=="unix") { ## bugfix for Wayland rstudio clients
  Sys.setenv(QT_QPA_PLATFORM="")
}
if(is.null(shinytest:::find_phantom())) stop("PhantomJS is not available")
