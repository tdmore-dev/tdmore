
#' Automatically plot a tdmorefit object
#'
#' The variable plotted on the y-axis is the first column with residual error that is available in the predicted result.
#' You can easily influence this by providing a `newdata` argument.
#'
#' @param x tdmorefit object
#' @param ... extra arguments
#' @importFrom ggplot2 autolayer
#' @param color color to use for the line
#' @param fill fill color to use for the ribbon
#' @param alpha alpha level for the ribbon
#' @param se.fit plot a confidence band around the prediction?
#' @param level level for the confidence band
#' @param na.rm should NA be removed
#' @export
#' @name autoplot.tdmore
autolayer.tdmorefit <- function(x, color="tomato1", fill=color, alpha=0.2, se.fit=T, level=0.95, ..., na.rm=TRUE) {
  ## Determine the mapping for these layers
  available <- stats::predict(x, newdata=numeric(), se.fit=se.fit, level=level) #these are all columns that are available
  newdata <- list(...)$newdata
  if(!is.null(newdata) && is.data.frame(newdata)) {
    available <- available[, colnames(newdata)]
  }
  mapping_line <- getAes(x, c("x","y"), data=available, as.aes=T)
  mapping_ribbon <- getAes(x, c("x", "ymin", "ymax"), data=available, as.aes=T)

  z <- list()

  z <- c(z, predictionLayer(mapping=mapping_line, data=x, geom="line", color=color, se.fit=se.fit, level=level, ..., na.rm=na.rm) ) #aes mapping determined automatically
  if(se.fit) z <- c(z, predictionLayer(mapping=mapping_ribbon, data=x, geom="ribbon", fill=fill, alpha=alpha, se.fit=se.fit, level=level, ..., na.rm=na.rm) )

  z <- c(z, autolayer_observed(x, mapping_line, available, na.rm))

  z
}

autolayer_observed <- function(x, mapping_line, available, na.rm) {
  yName <- getAes(x, "y", data=available, as.aes=F)$y
  observed <- model.frame(x)
  if(! yName %in% colnames(observed)) return(NULL)
  yValues <- observed[ , yName, drop=TRUE]

  if(
    (!na.rm && length(yValues) > 0)
    ||
    (na.rm && length( stats::na.omit(yValues) ) > 0 )
  ) {
    ## OK, we can also display observed points!
    return( geom_point(mapping=mapping_line, data=observed, na.rm=na.rm)) # add default na.rm=TRUE
  }
}

#' Automatically create the required layers for a recommendation object
#' @importFrom ggplot2 autolayer
#' @inheritParams autoplot.tdmore
#' @param x recommendation
#' @param ... extra arguments
#' @export
autolayer.recommendation <- function(x, color="green", fill=color, alpha=0.2, se.fit=T, level=0.95, ...) {
  tNames <- colnames(x$target)
  c(
    autolayer.tdmorefit(x$tdmorefit, regimen=x$regimen, color=color, fill=fill, alpha=alpha, se.fit=se.fit, level=level, ...),
    list(
      geom_point(mapping=aes_string(x=tNames[1], y=tNames[2]), data=x$target, shape=13, size=4)
    )
  )
}

#' Plot a tdmorefit_mixture object.
#'
#' @param x the tdmorefit object
#' @param population should I plot the population prediction (in blue) ?
#' @param fit should I plot the fit (in red) ?
#' @param ... passed on to `predict`
#'
#' @return a ggplot object with the fitted individual curve, the 95% CI of this curve, the population prediction (with between-subject variability) and the observed data
#'
#' @importFrom ggplot2 autoplot
#' @export
autoplot.tdmorefit_mixture <- function(x, population=TRUE, fit=TRUE, ...) {
  autoplot.tdmorefit(x, population=population, fit=fit, ...)
}

#' Plot a tdmorefit object.
#'
#' @param x the tdmorefit object
#' @param population should I plot the population prediction (in blue) ?
#' @param fit should I plot the fit (in red) ?
#' @param ... passed on to `predict`
#'
#' @return a ggplot object with the fitted individual curve, the 95% CI of this curve, the population prediction (with between-subject variability) and the observed data
#'
#' @importFrom ggplot2 autoplot
#' @export
autoplot.tdmorefit <- function(x, population=TRUE, fit=TRUE, ...) {
  z <- ggplot(data=x, ...)

  z <- z + ggplot2::scale_x_continuous() + ggplot2::scale_y_continuous()

  if(population) z <- z + autolayer( as.population(x), color="skyblue1", ... ) #population
  if(fit) z <- z + autolayer(x, ...)

  z
}

#' Plot a tdmore object as the population
#' @param x a tdmore object
#' @param regimen the regimen
#' @param covariates the covariates
#' @param newdata passed to `predict`
#' @param ... extra arguments
#'
#' @importFrom ggplot2 autoplot
#' @export
autoplot.tdmore <- function(x, regimen=NULL, covariates=NULL, newdata=NULL, ...) {
  x <- x %>%
    estimate(regimen = regimen, covariates = covariates)
  if(is.null(newdata)) newdata <- timeRange(x)
  z <- autoplot.tdmorefit(x, ..., newdata=newdata, population=TRUE, fit=FALSE)
  z
}

#' Determine a timerange suitable for plotting each aspect of a tdmorefit
#' @param ipred tdmorefit object to calculate the time range
#' @param regimen the regimen, or not specified to use the original regimen
#' @param observed observed values, or not specified to use the original observed values
#' @param N number of interpolation points
#'
#' @noRd
timeRange <- function(ipred=NULL, regimen=NULL, observed=NULL, N=101) {
  if(is.null(regimen) && !is.null(ipred)) regimen <- ipred$regimen
  if(is.null(observed) && !is.null(ipred)) observed <- ipred$observed
  maxTime <- computeMaxTime(regimen=regimen, observed=observed)
  times <- seq(0, maxTime, length.out=N)
  #ensure the regimen/ctrough are included in the points
  times <- unique(c(times, regimen$TIME, regimen$TIME[regimen$TIME > 0]-.Machine$double.eps, observed$TIME)) %>% sort
  times
}


#' Compute the maximum interesting time value based on the specified regimen and observed data.
#'
#' @param regimen the specified regimen
#' @param observed the observed dataframe
#'
#' @details the maximum of all times in observed, and regimen (+one additional Interdose-interval, should this be specified)
#'
#' @return the tmax value, numeric
#' @noRd
computeMaxTime <- function(regimen, observed=NULL) {
  values <- c(0)
  if (!is.null(regimen)) {
    values <- c(values, regimen$TIME)
    if(all(c("ADDL", "II") %in% names(regimen))) values <- c(values, regimen$TIME + (regimen$ADDL+1) * regimen$II)
  }
  if (!is.null(observed)) {
    values <- c(values, observed$TIME)
  }
  return(max(values, na.rm=TRUE))
}

#' @importFrom graphics plot
#' @export
plot.tdmorefit <- autoplot.tdmorefit


#' @importFrom graphics plot
#' @export
plot.tdmore <- autoplot.tdmore

#' @importFrom graphics plot
#' @export
plot.tdmorefit_mixture <- autoplot.tdmorefit_mixture

#' Generate a plot of the parameters
#'
#' @inheritParams predict.tdmorefit
#' @importFrom ggplot2 aes_
#' @param x tdmorefit object
#' @param vars the variables to plot, as a character vector. If missing, this is taken from the observed_variables metadata
#'
#' @export
parameterPlot.tdmorefit <- function(x, newdata, vars, ...) {
  z <- ggplot(data=x, ...)
  z <- z + ggplot2::scale_x_continuous() + ggplot2::scale_y_continuous()

  observedVariables <- getMetadataByClass(x$tdmore, "tdmore_observed_variables")
  if(missing(vars)) vars <- observedVariables$variables
  if(is.null(vars) || length(vars) == 0) stop("No observed_variables metadata specified!")
  df <- stats::predict(x, newdata=newdata) %>% tidyr::pivot_longer(cols=!!vars, values_to="IPRED")
  dfBase <- stats::predict(as.population(x), newdata=newdata) %>% tidyr::pivot_longer(cols=!!vars, values_to="PRED")
  df$PRED <- dfBase$PRED
  #df <- dplyr::bind_cols(df, dfBase)

  for(i in observedVariables$variables) {
    z <- z + geom_line(aes_(x=~TIME, y=~IPRED/PRED - 1, linetype=i),
                       data=dplyr::filter(df, .data$name == !!i) )
  }

  z + labs(x="Time" , y="Deviation from typical value ", linetype="")
}
