
plot.tdmorefit <- function(tdmorefit, newdata=NULL, .progress="none") {
  if(is.null(newdata)) {
    newdata <- data.frame(TIME=seq(0, max(tdmorefit$observed$TIME), length.out=100))
    for(i in colnames(tdmorefit$observed)) {
      if(i == "TIME") next
      newdata[, i] <- NA
    }
  }
  oNames <- names(newdata)
  oNames <- oNames[oNames != "TIME"]
  melt <- function(x, measure.vars=oNames, se=FALSE) {
    vars <- measure.vars
    if( se ) {
      for(i in c(".upper", ".lower")) vars <-c(vars, paste0(measure.vars, i))
    }
    tmp <- reshape::melt(x, id.vars="TIME", measure.vars=vars)
    if(se) {
      result <- subset(tmp, variable %in% measure.vars)
      result$value.upper <- subset(tmp, variable %in% paste0(measure.vars, ".upper"))$value
      result$value.lower <- subset(tmp, variable %in% paste0(measure.vars, ".lower"))$value
    } else {
      result <- tmp
    }
    result
  }

  ipred <- tdmorefit %>% predict(newdata) %>% melt
  ipredre <- tdmorefit %>% predict.tdmorefit(newdata, se.fit=TRUE, level=0.95, .progress=.progress) %>% melt(se=TRUE)
  pred <- estimate(tdmorefit$tdmore, regimen=tdmorefit$regimen) %>% predict(newdata) %>% melt
  obs <- model.frame.tdmorefit(tdmorefit, se=TRUE, level=0.95) %>% melt(se=TRUE)

  z <- ggplot(mapping=aes(x=TIME, y=value)) +
    geom_line(aes(color="Fit"), data=ipred) +
    geom_ribbon(aes(fill="Fit (95% CI)", ymin=value.lower, ymax=value.upper), data=ipredre, alpha=0.3)+
    geom_line(aes(color="Population"), data=pred) +
    geom_point(data=obs) +
    geom_errorbar(aes(ymax=value.upper, ymin=value.lower), data=obs)
  z
}
