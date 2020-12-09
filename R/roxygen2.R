#' The @@engine tag marks a function as part of the internal engine.
#' It is not meant for external use, and should therefore not generate
#' external-facing documentation or be included in the NAMESPACE
#'
#' @importFrom roxygen2 roxy_tag_parse
#' @export
#' @noRd
roxy_tag_parse.roxy_tag_engine <- function(x) {
  ## The noRd tag will stop Rd generation
  # the documentation will not be included in
  # the pkgdown site either.
  replacement <- roxygen2::roxy_tag_parse(
    roxygen2::roxy_tag(tag="noRd", raw="")
  )
  return(replacement)
}
