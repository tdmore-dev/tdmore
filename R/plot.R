
#' Plot a tdmorefit object.
#'
#' @param x the tdmorefit object
#' @param newdata a data.frame with at least TIME and any other columns to plot, NULL to plot all columns from the original observed data between time 0 and max(observationTime) or a numeric vector of times
#' @param se.fit add a curve for the confidence interval around the fit
#' @param mc.maxpts maximum number of points to use for the monte carlo fit
#' @param .progress either "none" or "text" to see calculation progress of monte carlo simulations
#' @param ... ignored
#'
#' @return a ggplot object with the fitted individual curve, the 95% CI of this curve, the population prediction (with between-subject variability) and the observed data
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon geom_point geom_errorbar labs
#' @importFrom stats predict
#' @export
plot.tdmorefit <- function(x, newdata=NULL, se.fit=TRUE, mc.maxpts=100, .progress="none", ...) {
  tdmorefit <- x
  newdata <- processNewData(newdata, tdmorefit)
  args <- list(...)
  populationArg <- unlist(args['population'])
  population <- if (length(populationArg)==0) F else populationArg

  ipred <- tdmorefit %>% predict(newdata) %>% meltPredictions()
  if (se.fit) {
    ipredre <- tdmorefit %>% predict.tdmorefit(newdata, se.fit=T, level=0.95, mc.maxpts = mc.maxpts, .progress=.progress) %>% meltPredictions(se=T)
  }
  yVars <- colnames(newdata)[colnames(newdata) != "TIME"]

  if (population) {
    # No need to compute PRED because IPRED = PRED
    plot <- ggplot(mapping=aes_string(x="TIME", y="value", group="variable")) + geom_line(color=blue(), data=ipred)
    if (se.fit) plot <- plot + geom_ribbon(fill=blue(), aes_string(ymin="value.lower", ymax="value.upper"), data=ipredre, alpha=0.05)

  } else {
    # Compute PRED
    pred <- estimate(tdmorefit$tdmore, regimen=tdmorefit$regimen, covariates=tdmorefit$covariates) %>% predict(newdata) %>% meltPredictions()
    if (se.fit) {
      predre <- estimate(tdmorefit$tdmore, regimen=tdmorefit$regimen, covariates=tdmorefit$covariates) %>% predict(newdata, se.fit=T) %>% meltPredictions(se=T)
    }

    obs <- model.frame.tdmorefit(tdmorefit) %>% meltPredictions()
    obs <- subset(obs, obs$variable %in% yVars & !is.na(obs$value))

    plot <- ggplot(mapping=aes_string(x="TIME", y="value", group="variable"))
    plot <- plot + geom_line(color=blue(), data=pred) + geom_point(data=obs)
    plot <- plot + geom_line(color=red(), data=ipred) + geom_point(data=obs)
    if (se.fit) plot <- plot + geom_ribbon(fill=blue(), aes_string(ymin="value.lower", ymax="value.upper"), data=predre, alpha=0.05)
    if (se.fit) plot <- plot + geom_ribbon(fill=red(), aes_string(ymin="value.lower", ymax="value.upper"), data=ipredre, alpha=0.04)
  }

  plot <- plot + labs(y=paste(yVars, collapse = " / "))
  return(plot)
}

#' Plot a tdmore object (typical profile of the population and between subject variability).
#'
#' @param x the tdmore object
#' @param regimen the regimen to be predicted
#' @param covariates covariates
#' @param newdata a data.frame with at least TIME and any other columns to plot, NULL to plot all columns from the original observed data between time 0 and max(observationTime) or a numeric vector of times
#' @param bsv add between subject variability, enabled by default
#' @param mc.subjects number of subjects in the monte carlo simulation
#' @param .progress either "none" or "text" to see calculation progress of monte carlo simulations
#' @param ... ignored
#'
#' @return a ggplot object with the the population prediction
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon geom_point geom_errorbar labs
#' @importFrom stats predict
#' @export
plot.tdmore <- function(x, regimen, covariates=NULL, newdata=NULL, bsv=TRUE, mc.subjects=100, .progress="none", ...) {
  tdmore <- x
  tdmorefit <- estimate(tdmore = tdmore, regimen = regimen, covariates = covariates)
  plot <- plot(tdmorefit, newdata=newdata, se.fit=bsv, mc.maxpts=mc.subjects, .progress=.progress, population=TRUE)
  return(plot)
}

processNewData <- function(newdata, tdmorefit) {
  if (is.null(newdata)) {
    newdata <- data.frame(TIME = seq(0, computeTmax(tdmorefit$regimen, tdmorefit$observed), length.out = 100))
    newdata <- addObservedVariables(newdata, tdmorefit$tdmore)
  }
  if (is.numeric(newdata)) {
    newdata <- data.frame(TIME = newdata)
    newdata <- addObservedVariables(newdata, tdmorefit$tdmore)
  }
  return(newdata)
}

addObservedVariables <- function(newdata, tdmore) {
  for (err in tdmore$res_var) {
    newdata[, err$var] <- NA
  }
  return(newdata)
}

blue <- function() {
  return("steelblue2")
}

red <- function() {
  return("tomato1")
}
