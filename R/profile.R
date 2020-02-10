
#' Plot a tdmoreprofile object.
#'
#' @param x the tdmoreprofile object
#' @param parameters the parameters to plot
#' @param contour use the 'geom_contour' method to draw the profile, logical value
#' @param raster use the 'raster method' to draw the profile, logical value
#' @param ... unused arguments
#' @importFrom ggplot2 ggplot aes_string geom_raster geom_contour geom_point geom_line labs
#' @export
plot.tdmoreprofile <- function(x, parameters=NULL, contour=T, raster=T, ...) {
  profile <- x
  tdmorefit <- profile$tdmorefit

  if(is.null(parameters)) {
    selectedParameters <- profile$profiledParameters
  } else {
    selectedParameters <- parameters
  }

  stopifnot(all(selectedParameters %in% colnames(profile$profile)))
  length <- length(selectedParameters)

  if(length > 2) {
    stop("Cannot produce plot in more than 2 dimensions... Use tdmore::shinyProfile() to explore the log-likelihood space instead")
  } else if(length == 2) {
    if(length(profile$profiledParameters) > 2) {
      warning(paste0("More than 2 profiled parameters, leading to multiple logLik-scores for a given parameter combination.\n",
           "Execute the profile() again with just 2 parameters, or use facets to ensure unique logLik scores per point."))
    }
    plot <- ggplot(data = profile$profile, aes_string(x=selectedParameters[1], y=selectedParameters[2], z="exp(logLik)"))
    if(raster) {plot <- plot + geom_raster(aes_string(fill = "exp(logLik)"))}
    if(contour){plot <- plot + geom_contour(colour="grey")}
    res <- coef(profile$tdmorefit)
    plot + geom_point(x=c(res[[selectedParameters[1]]]), y=c(res[[selectedParameters[2]]]), colour="white", inherit.aes=FALSE)

  } else if(length == 1) {
    notSelected <- setdiff(profile$profiledParameters, selectedParameters)
    group <- paste0("interaction(", paste(notSelected, collapse=","), ")")
    if(length(notSelected)==0) group <- NULL
    ggplot(data = profile$profile, aes_string(x=selectedParameters[1], y="exp(logLik)", group=group)) +
      geom_line()
  } else {
    stop("No plot can be produced. There should be at least 1 or 2 parameters to draw.")
  }
}

