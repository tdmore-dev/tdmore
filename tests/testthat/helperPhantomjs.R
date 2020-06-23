setup({
  path <- shinytest:::find_phantom()
  if(.Platform$OS.type=="unix") { ## bugfix for Wayland rstudio clients
    Sys.setenv(QT_QPA_PLATFORM="")
  }
  if(is.null(path)) stop("PhantomJS is not available")
})
