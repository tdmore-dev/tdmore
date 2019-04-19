#' Metadata is a generic function to append metadata to a TDMore model.
#'
#' @param ... extra arguments
#' @return an object of class tdmore
#' @export
metadata <- function(...) {
  UseMethod("metadata")
}

#' Default function, catching unsupported types of metadata
#'
#' @param x anything
#' @param ... ignored
metadata.default <- function(x, ...) {
  stop("Metadata of class ", class(x), " not supported.")
}

#' Build a covariate metadata.
#'
#' @param x a tdmore_covariate object
#' @param ... ignored
#' @return a metadata structure
#' @export
metadata.tdmore_covariate <- function(x, ...) {
  structure(list(
    covariate=x
  ), class="metadata")
}

#' Build an output metadata.
#'
#' @param x a tdmore_output object
#' @param ... ignored
#' @return a metadata structure
#' @export
metadata.tdmore_output <- function(x, ...) {
  structure(list(
    output=x
  ), class="metadata")
}

#' Append metadata to a TDMore model.
#'
#' @param tdmore a tdmore_output object
#' @param ... metadata
#' @return a tdmore model
#' @export
metadata.tdmore <- function(tdmore, ...) {
  args <- list(...)
  index <- length(tdmore$metadata)
  for(arg in args) {
    index <- index + 1
    tdmore$metadata[[index]] <- arg
  }
  return(tdmore)
}

#' Create a new covariate.
#'
#' @param name the covariate name
#' @param label the covariate label
#' @param unit the unit, can be null
#' @param min the minimum value, can be null
#' @param max the maximum value, can be null
#' @param choices a list of choices, e.g. list(Fast=0, Slow=1)
#' @return a tdmore_covariate object
#' @export
covariate <- function(name, label, unit=NULL, min=NULL, max=NULL, choices=NULL) {
  structure(list(
    name=name,
    label=label,
    unit=unit,
    min=min,
    max=max,
    choices=choices
  ), class="tdmore_covariate")
}

#' Create a new output
#'
#' @param name the output name
#' @param label the output label
#' @param unit the unit, can be null
#' @return a tdmore_output object
#' @export
output <- function(name, label, unit=NULL) {
  structure(list(
    name=name,
    label=label,
    unit=unit
  ), class="tdmore_output")
}
