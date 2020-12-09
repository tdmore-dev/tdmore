#' Get the dosing interval for the given formulation
#'
#' @param x formulation (character vector)
#' @param model a tdmore model with formulation metadata
#'
#' @examples
#' if(FALSE) {
#' regimen <- data.frame(TIME=0, AMT=5, FORM="CompoundA")
#' regimen$II <- getDosingInterval(regimen$FORM, model)
#' }
#' @export
#' @family metadata
getDosingInterval <- function(x, model) {
  if(length(x) == 0) return(numeric(0))
  if(length(x)>1) return( purrr::map_dbl(x, getDosingInterval, model) )
  form <- tdmore::getMetadataByName(model, x)
  if(is.null(form)) stop("Formulation `", x, "' is not defined in the model metadata")
  form$dosing_interval
}

#' Append metadata to a TDMore model
#'
#' @param tdmore a tdmore_output object
#' @param ... metadata
#' @return a tdmore model
#' @export
#' @family metadata
metadata <- function(tdmore, ...) {
  if(!is.tdmore(tdmore)) stop("Can only add metadata to an existing tdmore object")
  args <- list(...)
  index <- length(tdmore$metadata)
  for(arg in args) {
    index <- index + 1
    tdmore$metadata[[index]] <- arg
  }
  names(tdmore$metadata) <- vapply(tdmore$metadata, function(x){x$name}, FUN.VALUE=character(1))

  return(tdmore)
}

#' Create a new covariate object.
#'
#' @param name the covariate name
#' @param label the covariate label
#' @param unit the unit, can be null
#' @param min the minimum value, can be null
#' @param max the maximum value, can be null
#' @param choices a numeric named vector of choices, e.g. c(Fast=0, Slow=1)
#' @return a tdmore_covariate object
#' @export
#' @family metadata
covariate <- function(name, label, unit=NULL, min=NULL, max=NULL, choices=NULL) {
  if(!is.null(choices)) stopifnot( is.numeric(choices) )
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
#' @family metadata
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
#' @family metadata
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
#' @family metadata
target <- function(min, max) {
  structure(list(
    name="TARGET",
    min=min,
    max=max
  ), class="tdmore_target")
}

#' Create an 'observed variables' object.
#' This metadata indicates that these variables are interesting outputs of the model,
#' and should be included in a plot.
#'
#' @param variables variables to be observed
#' @return a tdmore_observed_variables object
#' @export
observed_variables <- function(variables) {
  stopifnot(is.character(variables))
  structure(list(
    name="OBSERVED_VARIABLES",
    variables=variables
  ), class="tdmore_observed_variables")
}

#' Get the observed variables from the tdmore model.
#'
#' @param model tdmore model
#' @return a character vector of observed variables
#' @export
getObservedVariables <- function(model) {
  metadata <- getMetadataByName(model, "OBSERVED_VARIABLES")
  if (is.null(metadata)) {
    return(c())
  } else {
    return(metadata$variables)
  }
}


#' @export
print.tdmore_output <- function(x, ...) {
  cat(paste0(x$label, " (", x$unit, ")"))
}

#' @export
print.tdmore_covariate <- function(x, ...) {
  cat(paste0(x$label, " (", x$unit, ")"))
}

#' @export
print.tdmore_formulation <- function(x, ...) {
  cat(x$unit)
}

#' @export
print.tdmore_target <- function(x, ...) {
  cat(paste0("Target: [", x$min, ",", x$max, "]"))
}

#' Get metadata (output, formulation or covariate) by name.
#'
#' @param tdmore the tdmore model
#' @param metaName the name of the metadata
#' @param all whether to return only a single element, or all of them (in a list)
#' @return the first metadata item that matches the given name
#' @export
#' @family metadata
getMetadataByName <- function(tdmore, metaName, all=FALSE) {
  if(is.null(metaName)) return(NULL)
  hasMetadata <- function(x) {if ('name' %in% names(x)) x$name==metaName else FALSE}
  if(length(tdmore$metadata)==0) return(NULL) #no metadata
  i <- vapply(tdmore$metadata, hasMetadata, FUN.VALUE=logical(1))
  if(any(i)) {
    if(all) tdmore$metadata[i] else tdmore$metadata[[which.max(i)]]
  } else {
    NULL
  }
}

#' Get metadata by class.
#'
#' @param tdmore the tdmore model
#' @param metaClass the class of the metadata
#' @return the metadata items that matches the given class
#' @param all whether to return only a single element, or all of them (in a list)
#' @export
#' @family metadata
getMetadataByClass <- function(tdmore, metaClass, all=FALSE) {
  if(is.null(metaClass)) return(NULL)
  if(length(tdmore$metadata)==0) return(NULL) #no metadata
  hasMetadata <- function(x) {inherits(x, metaClass) }
  i <- vapply(tdmore$metadata, hasMetadata, FUN.VALUE=logical(1))
  if(any(i)) {
    if(all) tdmore$metadata[i] else tdmore$metadata[[which.max(i)]]
  } else {
    NULL
  }
}
