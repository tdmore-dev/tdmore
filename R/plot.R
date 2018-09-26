
#' Plot a tdmorefit object
#'
#' @param x the tdmorefit object
#' @param newdata a data.frame with at least TIME and any other columns to plot,
#' or NULL to plot all columns from the original observed data between
#' time 0 and max(observationTime),
#' or a numeric vector of times
#' @param se.fit Add a curve for the confidence interval around the fit
#' @param se.obs Add error bars for the residual error around the observations
#' @param mc.maxpts Maximum number of points to use for the monte carlo fit
#' @param .progress either "none" or "text" to see calculation progress of monte carlo simulations
#' @param ... ignored
#'
#' @return a ggplot object with the fitted individual curve, the 95% ci of this curve, the population prediction, the observed data and the residual error around the observed data
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon geom_point geom_errorbar
#' @importFrom stats predict
#' @export
plot.tdmorefit <- function(x, newdata=NULL, se.fit=TRUE, se.obs=TRUE, mc.maxpts=100, .progress="none", ...) {
  tdmorefit <- x
  if(is.null(newdata)) {
    tmax <- max(0, tdmorefit$observed$TIME,
                tdmorefit$regimen$TIME,
                tdmorefit$regimen$TIME + tdmorefit$regimen$ADDL * tdmorefit$regimen$II,
                na.rm=TRUE)
    newdata <- data.frame(TIME=seq(0, tmax, length.out=100))
    for(i in colnames(tdmorefit$observed)) {
      if(i == "TIME") next
      newdata[, i] <- NA
    }
  }
  if(is.numeric(newdata)) {
    #TODO
  }
  #oNames <- names(newdata)
  #oNames <- oNames[oNames != "TIME"]
  melt <- function(x, se=FALSE) {
    measure.vars <- colnames(x)
    measure.vars <- measure.vars[measure.vars != "TIME"]
    vars <- measure.vars
    if( se ) {
      for(i in c(".upper", ".lower")) vars <-c(vars, paste0(measure.vars, i))
    }
    tmp <- reshape::melt(x, id.vars="TIME")
    if(se) {
      result <- tmp[ tmp$variable %in% measure.vars , ]
      result$value.upper <- tmp$value[ tmp$variable %in% paste0(measure.vars, ".upper") ]
      result$value.lower <- tmp$value[ tmp$variable %in% paste0(measure.vars, ".lower") ]
    } else {
      result <- tmp
    }
    result
  }

  ipred <- tdmorefit %>% predict(newdata) %>% melt
  ipredre <- tdmorefit %>% predict.tdmorefit(newdata, se.fit=TRUE, level=0.95, mc.maxpts = mc.maxpts, .progress=.progress) %>% melt(se=TRUE)
  pred <- estimate(tdmorefit$tdmore, regimen=tdmorefit$regimen) %>% predict(newdata) %>% melt
  obs <- model.frame.tdmorefit(tdmorefit, se=TRUE, level=0.95) %>% melt(se=TRUE)

  z <- ggplot(mapping=aes_string(x="TIME", y="value")) +
    geom_line(color=1, data=ipred)
  if(se.fit) z <- z + geom_ribbon(fill=1, aes_string(ymin="value.lower", ymax="value.upper"), data=ipredre, alpha=0.3)
  z <- z +
    geom_line(color=2, data=pred) +
    geom_point(data=obs)
  if(se.obs) z <- z + geom_errorbar(aes_string(ymax="value.upper", ymin="value.lower"), data=obs)
  z
}

