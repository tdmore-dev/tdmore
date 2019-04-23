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

#' Build a dose metadata.
#'
#' @param x a tdmore_dose object
#' @param ... ignored
#' @return a metadata structure
#' @export
metadata.tdmore_dose <- function(x, ...) {
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

#' Create a new covariate object.
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

#' Create a new output object.
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

#' Create a new dose object.
#'
#'
#' @param unit the unit
#' @param name by default 'DOSE'
#' @return a tdmore_dose object
#' @export
dose <- function(unit, name="DOSE") {
  structure(list(
    name=name,
    unit=unit
  ), class="tdmore_dose")
}

#' Tdmore_output to string method.
#'
#' @param x a tdmore object
#' @param ... ignored
#'
#' @export
toString.tdmore_output <- function(x, ...) {
  return(paste0(x$label, " (", x$unit, ")"))
}

#' Tdmore_covariate to string method.
#'
#' @param x a tdmore object
#' @param ... ignored
#'
#' @export
toString.tdmore_covariate <- function(x, ...) {
  return(paste0(x$label, " (", x$unit, ")"))
}

#' Tdmore_dose to string method.
#'
#' @param x a tdmore object
#' @param ... ignored
#'
#' @export
toString.tdmore_dose <- function(x, ...) {
  return(x$unit)
}

#' Get metadata (output or covariate) by name.
#'
#' @param tdmore the tdmore model
#' @param metaName the name of the metadata
#' @return the first metadata item that matches the given name
#' @export
getMetadataByName <- function(tdmore, metaName) {
  hasMetadata <- function(x) {if ('name' %in% names(x)) x$name==metaName else FALSE}
  results <- tdmore$metadata[sapply(tdmore$metadata, hasMetadata)]
  if(length(results) > 0) {
    return(results[[1]])
  } else {
    return(NULL)
  }
}
