ggplot.lm <- function(model, ...) {
  base <- ggplot(fortify(model), ...)
  base$model <- model
  class(base) <- c("gglm", class(base))
  base
}

stat_predict <- function(data) { structure(list(data=data), class="StatPredict") }

ggplot_add.StatPredict <- function(object, plot, object_name) {
  if(!inherits(plot, "gglm")) stop("Can only add prediction to gglm objects")
  myData <- object$data
  myData$prediction <- predict(plot$model, newdata=object$data)
  plot + geom_line(aes(y=prediction), data=as.data.frame(myData))
}


mod <- lm(mpg ~ wt + hp, data = mtcars)
ggplot(mod, aes(x=wt, y=mpg)) +
  geom_point(aes(size=hp)) +
  stat_predict(data=data.frame(wt=seq(1, 6, length.out=100), hp=50)) +
  stat_predict(data=data.frame(wt=seq(1, 6, length.out=100), hp=300))
