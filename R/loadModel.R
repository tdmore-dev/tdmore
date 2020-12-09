#' This searches for the given model name in a directory
#' @param name model name (filename without `.r`)
#' @param dir model directory; uses the package `models` directory by default
#' @param pkg the package in which to search for a `models` directory
#' @return a tdmore model object
#' @export
getModel <- function(name=defaultModel(), dir=system.file(package=pkg, "models"), pkg="tdmore") {
  file <- dir(path=dir, pattern=paste0(name, "\\.[rR]$"))
  if(length(file)==0) stop("Model `", name, "' not found in directory `", dir, "'")
  env <- new.env()
  result <- source(file.path(dir, file), local=env, keep.source=TRUE)
  result$value
}

#' Get the first available model name in the given directory
#' @inheritParams getModel
#' @return character string with model filename
#' @export
defaultModel <- function(dir=system.file(package=pkg, "models"), pkg="tdmore"){
  available <- listModels(dir=dir)
  if(length(available)==0) stop("Not a single model available in ", dir)
  available[1]
}



#' This lists the models in a directory
#' @inheritParams getModel
#' @return a character vector of available model names
#' @noRd
listModels <- function(dir=system.file(package=pkg, "models"), pkg="tdmore") {
  available <- dir(path=dir, pattern="\\.[rR]$")
  #if(length(available)==0) stop("No models available in directory ",dir)
  sub("\\.[rR]$", "", available) #remove last part
}
