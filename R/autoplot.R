#' Plot a tdmorefit object.
#'
#' @param x the tdmorefit object
#' @param ... ignored
#'
#' @return a ggplot object with the fitted individual curve, the 95% ci of this curve, the population prediction, the observed data and the residual error around the observed data
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 ggplot aes aes_string geom_line geom_ribbon geom_point geom_errorbar labs
#' @importFrom stats predict
#' @export
autoplot.tdmorefit <- function(x, ...) {
  pop <- estimate(x$tdmore, regimen=x$regimen, covariates=x$covariates)
  tMax <- max(x$regimen$TIME)

  data=x$observed
  if(nrow(x$observed) == 0) data=model.frame(x$tdmore, data=NA) #simple data.frame with 1 row

  z1 <- ggplot(data=x, ..., mapping=aes_string(x="TIME"))
  if(nrow(x$observed) > 0) {
    firstVar <- as.name( colnames(x$observed)[2] )
    z1 <- z1 + geom_point(aes_(color="Observed", y=firstVar))
  }
  z1 +
    stat_predict(aes(color="A posteriori"), xlim=c(0, tMax), data=data) +
    stat_predict(aes(color="A priori"), tdmorefit=pop, xlim=c(0, tMax), data=data) +
    labs(x="Time", y="Value", color="Prediction", fill="Prediction")
}

#' Plot a tdmore object (typical profile (population) + between subject variability).
#'
#' @param x the tdmore object
#' @param regimen the regimen to be predicted
#' @param covariates covariates
#' @param ... ignored
#'
#' @return a ggplot object with the fitted individual curve, the 95% ci of this curve, the population prediction, the observed data and the residual error around the observed data
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon geom_point geom_errorbar labs autoplot
#' @importFrom stats predict
#' @export
autoplot.tdmore <- function(x, regimen=NULL, covariates=NULL, ...) {
  x <- estimate(x, regimen = regimen, covariates = covariates)
  tMax <- max(regimen$TIME)
  ggplot(data=x, mapping=aes_string(x="TIME")) +
    stat_predict(aes(color="A priori"), xlim=c(0, tMax), data=data.frame(TIME=0:tMax)) +
    labs(x="Time", y="Value", color="Prediction", fill="Prediction")
}
