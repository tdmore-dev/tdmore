#' Create an 'observed variables' object.
#' These variables will be predicted and shown in the panel below the timeline.
#' TODO: discuss if we keep this info as part of the UI metadata of the model.
#' Should we need a config tab to configure the UI metadata?
#' Do we want to persist this information in DB? Per model? Per patient? etc.
#'
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

#' Build an 'observed variables' metadata.
#'
#' @param x a tdmore_observed_variables object
#' @param ... ignored
#' @return a metadata structure
#' @export
metadata.tdmore_observed_variables <- function(x, ...) {
  structure(list(
    value=x
  ), class="metadata")
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
