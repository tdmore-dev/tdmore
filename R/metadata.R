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
    value=x
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
    value=x
  ), class="metadata")
}

#' Build a formulation metadata.
#'
#' @param x a tdmore_formulation object
#' @param ... ignored
#' @return a metadata structure
#' @export
metadata.tdmore_formulation <- function(x, ...) {
  structure(list(
    value=x
  ), class="metadata")
}

#' Build a target metadata.
#'
#' @param x a tdmore_target object
#' @param ... ignored
#' @return a metadata structure
#' @export
metadata.tdmore_target <- function(x, ...) {
  structure(list(
    value=x
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
#' @param default_value the default dose value, same unit as the one provided
#' @return a tdmore_output object
#' @export
output <- function(name, label, unit=NULL, default_value=1) {
  structure(list(
    name=name,
    label=label,
    unit=unit,
    default_value=default_value
  ), class="tdmore_output")
}

#' Create a new formulation object.
#'
#' @param name the formulation name
#' @param unit the unit
#' @param dosing_interval the dosing interval in hours, 24h is the default value
#' @param default_value the default dose value, same unit as the one provided
#' @param round_function the rounding function to round the recommended dose, by default doses are not rounded
#' @return a tdmore_formulation object
#' @export
formulation <- function(name, unit, dosing_interval=24, default_value=1, round_function=function(dose){dose}) {
  structure(list(
    name=name,
    unit=unit,
    dosing_interval=dosing_interval,
    default_value=default_value,
    round_function=round_function
  ), class="tdmore_formulation")
}

#' Create a new target.
#'
#' @param min the min value of the target, same unit as the one defined in output
#' @param max the max value of the target, same unit as the one defined in output
#' @return a tdmore_target object
#' @export
target <- function(min, max) {
  structure(list(
    name="TARGET",
    min=min,
    max=max
  ), class="tdmore_target")
}

#' Tdmore_output to string method.
#'
#' @param x a tdmore_output object
#' @param ... ignored
#'
#' @export
toString.tdmore_output <- function(x, ...) {
  return(paste0(x$label, " (", x$unit, ")"))
}

#' Tdmore_covariate to string method.
#'
#' @param x a tdmore_covariate object
#' @param ... ignored
#'
#' @export
toString.tdmore_covariate <- function(x, ...) {
  return(paste0(x$label, " (", x$unit, ")"))
}

#' Tdmore_formulation to string method.
#'
#' @param x a tdmore_dose object
#' @param ... ignored
#'
#' @export
toString.tdmore_formulation <- function(x, ...) {
  return(x$unit)
}

#' Tdmore_target to string method.
#'
#' @param x a tdmore_target object
#' @param ... ignored
#'
#' @export
toString.tdmore_target <- function(x, ...) {
  return(paste0("Target: [", x$min, ",", x$max, "]"))
}

#' Get metadata (output, formulation or covariate) by name.
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

#' Get metadata by class.
#'
#' @param tdmore the tdmore model
#' @param metaClass the class of the metadata
#' @return the metadata items that matches the given class
#' @export
getMetadataByClass <- function(tdmore, metaClass) {
  hasMetadata <- function(x) {inherits(x, metaClass) }
  results <- tdmore$metadata[sapply(tdmore$metadata, hasMetadata)]
  return(results)
}
