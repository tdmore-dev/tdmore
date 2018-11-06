
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
#' @importFrom graphics plot
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
#' @importFrom graphics plot
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


## TODO: Simplify all of the above; delegate to the code below.



#' Plot a tdmorefit object.
#'
#' @param x the tdmorefit object
#' @param newdata a data.frame with at least TIME and any other columns to plot, NULL to plot all columns from the original observed data between time 0 and max(observationTime) or a numeric vector of times
#' @param se.fit add a curve for the confidence interval around the fit
#' @param mc.maxpts maximum number of points to use for the monte carlo fit
#' @param .progress either "none" or "text" to see calculation progress of monte carlo simulations
#' @param se add residual error bounds, currently not compatible with se.fit
#' @param ... ignored
#'
#' @return a ggplot object with the fitted individual curve, the 95% ci of this curve, the population prediction, the observed data and the residual error around the observed data
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon geom_point geom_errorbar labs
#' @importFrom stats predict
#' @export
autoplot.tdmorefit <- function(x, newdata=NULL, ...) {
  ggplot(data=x, newdata=newdata) +
    pred_se(aes(fill="A priori")) +
    ipred_se(aes(fill="A posteriori")) +
    pred(aes(color="A priori")) +
    ipred(aes(fill="A posteriori")) +
    labs(x="Time", y="Value", color="Prediction", fill="Prediction")
}

#' Plot a tdmore object (typical profile (population) + between subject variability).
#'
#' @param x the tdmore object
#' @param regimen the regimen to be predicted
#' @param covariates covariates
#' @param newdata a data.frame with at least TIME and any other columns to plot, NULL to plot all columns from the original observed data between time 0 and max(observationTime) or a numeric vector of times
#' @param bsv add between subject variability
#' @param mc.subjects number of subjects in the monte carlo simulation
#' @param .progress either "none" or "text" to see calculation progress of monte carlo simulations
#' @param se add residual error bounds, currently not compatible with bsv
#' @param ... ignored
#'
#' @return a ggplot object with the fitted individual curve, the 95% ci of this curve, the population prediction, the observed data and the residual error around the observed data
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon geom_point geom_errorbar labs
#' @importFrom stats predict
#' @export
autoplot.tdmore <- function(x, regimen, covariates=NULL, newdata=NULL, ...) {
  tdmore <- x
  tdmorefit <- estimate(tdmore = tdmore, regimen = regimen, covariates = covariates)
  ggplot(data=tdmorefit, newdata=newdata) +
    pred_se(aes(fill="A priori")) +
    pred(aes(color="A priori")) +
    labs(x="Time", y="Value", color="Prediction", fill="Prediction")
}

## Design considerations
## Author: Ruben Faelens
## There are multiple ways to create a plotting library for tdmore
## 1) Create one single "god" plot functions, with a huge list of arguments
## 2) Split this up into several smaller functions. This is how xpose does it.
## Both of the above will always be annoying, because you need to learn a lot of extra helper functions...
##
## 3) Create helper functions that create simple ggplot layers. Disadvantage is that you must
## repeat arguments when you want to create multiple layers.
## 4) Create a separate S3 class for tdmoreplot objects. This is a ggplot object, but also a tdmoreplot object!
## Override the "+" to consistently adapt the created object in the right way.
## 5) Use the default extension mechanisms from ggplot (fortify, autoplot) in combination with the above
##
## So should we create a geom_ or a stat_ object?
## Basically, it does not really matter. See ?geom_histogram
## A stat_ object makes more sense, as it allows you to pick the geom you want.
## We will take our inspiration from stat_identity
##
##
## What does NOT work:
## 1) Use fortify to add extra attributes to the data.frame. The data.frame is already broken up when we are facetting,
## so any Geom or Stat* object would not see the original data.frame.
## Therefore, it is clear we need to calculate the new data.frame in ggplot_add. Only unfortunate thing:
## we cannot substitute data later on.
## However, as long as the data was created during ggplot_add, later facetting can happily split up that data.frame further!

#' @export
ggplot.tdmorefit <- function(data=NULL, mapping=aes(),
                             newdata=NULL, regimen=NULL, parameters=NULL, covariates=NULL,
                             ...,
                             environment=parent.frame()) {
  ggdata <- predict(data, newdata=newdata, regimen=regimen, parameters=parameters, covariates=covariates)
  base <- ggplot2::ggplot(ggdata, mapping, ..., environment)
  class(base) <- c("ggtdmore", class(base)) #needs to be the first one
  base$tdmore <- list(
    tdmorefit=data, newdata=newdata, regimen=regimen, parameters=parameters, covariates=covariates
  )
  base
}

## TODO: add ggplot functions for the other tdmore objects
## Ideally, we simply translate to the ggplot_tdmorefit structure.


## TODO: add the following functions
## ipred_se(geom="ribbon", fill="red")
## pred(geom="line", color="blue")
## pred_se(geom="ribbon", fill="blue", alpha=0.1)
## observed(geom="point")
## geom_target (eventually)
##

#' Create a population prediction.
#' @export
ipred <- function(mapping = NULL, data = NULL,
                  geom = "line", position = "identity",
                  color="steelblue2",
                  ...,
                  newdata=NULL,
                  regimen=NULL,
                  parameters=NULL,
                  covariates=NULL,
                  show.legend = NA,
                  inherit.aes = TRUE) {
  if(!is.null(data) && !is.tdmorefit(data)) stop("Only a tdmorefit object can be used as `data` for ipred")
  structure(
    list(
      data = list(
        tdmorefit=NULL,
        newdata=NULL,
        regimen=NULL,
        parameters=NULL,
        covariates=NULL
      ),
      mapping = mapping,
      stat = StatIPred,
      geom = geom,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        ...
      )
    ), class="IPred"
  )
}

#' @export
ggplot_add.IPred <- function(object, plot, object_name) {
  if(!is.ggtdmore(plot)) stop(object_name, " can only be used with tdmorefit ggplot objects. Did you remember to pass a tdmorefit object to ggplot() ?")

  tdmorefit <-  object$data$tdmorefit %||% plot$tdmore$tdmorefit
  newdata <-    object$data$newdata %||% plot$tdmore$newdata
  regimen <-    object$data$regimen %||% plot$tdmore$regimen
  parameters <- object$data$parameters %||% plot$tdmore$parameters
  covariates <- object$data$covariates %||% plot$tdmore$covariates

  object$data <- predict(tdmorefit, newdata=newdata, regimen=regimen, parameters=parameters, covariates=covariates)
  ## TODO: guess the y aesthetic
  if(is.null(object$mapping$y)) stop("No y aesthetic supplied; TODO: guess the y aesthetic from the error model")

  with(object, {
    plot + layer(
      data = data,
      mapping = mapping,
      stat = StatIPred,
      geom = geom,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = params
    )
  })
}

## How can we guess the y variable as well?
## We need the model for that ?
StatIPred <- ggplot2::ggproto("StatIPred", ggplot2::StatIdentity,
                              default_aes=aes(x=TIME),
                              required_aes = c("y")
)

is.ggtdmore <- function(a) {inherits(a, "ggtdmore")}
"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}


## TODO: Get inspired by the code below when creating ipred_se and pred_se
## Afterwards, you can remove it.
#' @importFrom ggplot2 ggplot_add aes
#' @export
ggplot_add.GeomPred <- function(object, plot, object_name) {
  if(!is.ggtdmore(plot)) stop("Can only add GeomPred to a ggtdmore object. Did you call ggplot() with a tdmore object?")

  ## If these arguments are missing, use the defaults from the original ggplot() call
  tdmorefit <- object$tdmorefit %||% plot$tdmore$tdmorefit
  newdata <- object$newdata %||% plot$tdmore$newdata
  regimen <- object$regimen %||% plot$tdmore$regimen
  parameters <- object$parameters %||% plot$tdmore$parameters
  covariates <- object$covariates %||% plot$tdmore$covariates

  ## First add the confidence interval
  if(object$se.fit) {
    data <-predict(tdmorefit, newdata=newdata, regimen=regimen, parameters=parameters, covariates=covariates,
                   se.fit=T, level=object$level)
    # TODO: guess the output, or use tidybayes

    mapping <- object$mapping %||% aes(x=TIME, ymin=CONC.lower, ymax=CONC.upper)
    # TODO: how do we work with color if it is also specified as an aesthetic?
    # TODO: What happens if
    args <- list(data=data, mapping=mapping, fill="red", alpha=0.2)
    args <- c(object$ribbonArgs, args)
    args <- args[!duplicated(names(args))] #drop duplicates
    layer <- do.call(geom_ribbon, args)
    plot <- plot + layer
  }

  data <- predict(tdmorefit, newdata=newdata, regimen=regimen, parameters=parameters, covariates=covariates)
  mapping <- object$mapping  %||% aes(x=TIME, y=CONC)# TODO: guess the output, or use tidybayes
  args <- list(data=data, mapping=mapping, color="red")
  args <- c(object$lineArgs, args)
  args <- args[!duplicated(names(args))] #drop duplicates
  layer <- do.call(geom_line, args)
  plot <- plot + layer

  plot
}
