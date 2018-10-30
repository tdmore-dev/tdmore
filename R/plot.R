
#' Plot a tdmorefit object.
#'
#' @param x the tdmorefit object
#' @param newdata a data.frame with at least TIME and any other columns to plot, NULL to plot all columns from the original observed data between time 0 and max(observationTime) or a numeric vector of times
#' @param vars additional variables to be plotted
#' @param se.fit add a curve for the confidence interval around the fit
#' @param mc.maxpts maximum number of points to use for the monte carlo fit
#' @param .progress either "none" or "text" to see calculation progress of monte carlo simulations
#' @param ... ignored
#'
#' @return a ggplot object with the fitted individual curve, the 95% ci of this curve, the population prediction, the observed data and the residual error around the observed data
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon geom_point geom_errorbar labs
#' @importFrom stats predict
#' @export
plot.tdmorefit <- function(x, newdata=NULL, vars=NULL, se.fit=TRUE, mc.maxpts=100, .progress="none", ...) {
  tdmorefit <- x
  newdata <- processNewData(newdata, regimen=tdmorefit$regimen, observed=tdmorefit$observed, vars=NULL)

  ipred <- tdmorefit %>% predict(newdata) %>% meltPredictions()
  if(se.fit) {
    ipredre <- tdmorefit %>% predict.tdmorefit(newdata, se.fit=T, level=0.95, mc.maxpts = mc.maxpts, .progress=.progress) %>% meltPredictions(se=T)
  }

  pred <- estimate(tdmorefit$tdmore, regimen=tdmorefit$regimen, covariates=tdmorefit$covariates) %>% predict(newdata) %>% meltPredictions()
  obs <- model.frame.tdmorefit(tdmorefit) %>% meltPredictions()

  plot <- ggplot(mapping=aes_string(x="TIME", y="value", group="variable")) + geom_line(color=red(), data=ipred)
  if(se.fit) plot <- plot + geom_ribbon(fill=red(), aes_string(ymin="value.lower", ymax="value.upper"), data=ipredre, alpha=0.03)
  plot <- plot + geom_line(color=blue(), data=pred) + geom_point(data=obs)
  yVars <- colnames(newdata)[colnames(newdata) != "TIME"]
  plot <- plot + labs(y=paste(yVars, collapse = " / "))
  return(plot)
}

#' Plot a tdmore object.
#'
#' @param x the tdmorefit object
#' @param regimen the regimen to be predicted
#' @param covariates covariates
#' @param newdata a data.frame with at least TIME and any other columns to plot, NULL to plot all columns from the original observed data between time 0 and max(observationTime) or a numeric vector of times
#' @param vars additional variables to be plotted
#' @param se add a curve for the confidence interval around the population prediction
#' @param mc.maxpts maximum number of points to use for the monte carlo fit
#' @param .progress either "none" or "text" to see calculation progress of monte carlo simulations
#' @param ... ignored
#'
#' @return a ggplot object with the fitted individual curve, the 95% ci of this curve, the population prediction, the observed data and the residual error around the observed data
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon geom_point geom_errorbar labs
#' @importFrom stats predict
#' @export
plot.tdmore <- function(x, regimen, covariates=NULL, newdata=NULL, vars=NULL, se=TRUE, mc.maxpts=100, .progress="none", ...) {
  tdmore <- x
  newdata <- processNewData(newdata, regimen, observed=NULL, vars)
  predict <- predict.tdmore(tdmore, newdata, regimen=regimen, covariates=covariates)
  pred <- tdmore %>% predict.tdmore(newdata, regimen=regimen, covariates=covariates) %>% meltPredictions()
  if(se) {
    predre <- tdmore %>% predict.tdmore(newdata, regimen=regimen, covariates=covariates, se=T, level=0.95, mc.maxpts = mc.maxpts, .progress=.progress) %>% meltPredictions(se=T)
  }

  plot <- ggplot(mapping=aes_string(x="TIME", y="value", group="variable")) + geom_line(color=blue(), data=pred)
  if(se) plot <- plot + geom_ribbon(fill=blue(), aes_string(ymin="value.lower", ymax="value.upper"), data=predre, alpha=0.06)
  yVars <- colnames(newdata)[colnames(newdata) != "TIME"]
  plot <- plot + labs(y=paste(yVars, collapse = " / "))
  return(plot)
}

processNewData <- function(newdata, regimen, observed, vars) {
  if (is.null(newdata)) {
    newdata <- data.frame(TIME = seq(0, computeTmax(regimen, observed), length.out = 100))
    newdata <- addObservedVariables(newdata, observed, vars)
  }
  if (is.numeric(newdata)) {
    newdata <- data.frame(TIME = newdata)
    newdata <- addObservedVariables(newdata, observed, vars)
  }
  return(newdata)
}

addObservedVariables <- function(newdata, observed, vars) {
  observedVariables <- c()
  if (!is.null(observed)) {
    observedVariables <- c(observedVariables, colnames(observed))
  }
  if (!is.null(vars)) {
    observedVariables <- c(observedVariables, vars)
  }

  for (i in unique(observedVariables)) {
    if (i == "TIME")
      next
    newdata[, i] <- NA
  }
  return(newdata)
}

blue <- function() {
  return("steelblue2")
}

red <- function() {
  return("tomato1")
}
