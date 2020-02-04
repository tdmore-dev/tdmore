
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

  if(length >= 2) {
    plot <- ggplot(data = profile$profile, aes_string(x=selectedParameters[1], y=selectedParameters[2], z="exp(logLik)"))
    if(raster) {plot <- plot + geom_raster(aes_string(fill = "exp(logLik)"))}
    if(contour){plot <- plot + geom_contour(colour="grey")}
    res <- coef(profile$tdmorefit)
    plot + geom_point(x=c(res[[selectedParameters[1]]]), y=c(res[[selectedParameters[2]]]), colour="white", inherit.aes=FALSE)

  } else if(length == 1) {
    ggplot(data = profile$profile, aes_string(x=selectedParameters[1], y="exp(logLik)")) +
      geom_line()
  } else {
    warning("No plot can be produced. There should be at least 1 or 2 parameters to draw.")
  }
}

