setup({
  path <- shinytest:::find_phantom()
  if(is.null(path)) stop("PhantomJS is not available")
})
