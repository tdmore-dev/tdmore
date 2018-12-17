## Test
rm(list=ls(all=TRUE))

fortify.lm2 <- function(model, data, ...) {
  base <- NextMethod()
  attr(base, 'model') <- model
  class(base) <- c("lm.data.frame", class(base))
  base
}

stat_predict <- function(data) {
  layer(geom="line", stat="identity", position="identity", data=data, layer_class=LMPredictLayer)
}
LMPredictLayer <- ggproto(
  "LMPredictLayer",
  ggplot2:::Layer,
  setup_layer = function(self, data, plot) {
    model <- attr(plot$data, 'model')
    terms <- all.vars(formula(model))
    data[, terms[1]] <- predict(model, newdata=data)
    data
  }
)

mod <- lm(mpg ~ wt + hp, data = mtcars)
class(mod) <- c("lm2", class(mod)) # to make sure our fortify method gets called first
z1 <- ggplot(mod, aes(x=wt, y=mpg)) +
  geom_point(aes(size=hp)) +
  stat_predict(data=data.frame(wt=seq(1, 6, length.out=100), hp=50)) +
  stat_predict(data=data.frame(wt=seq(1, 6, length.out=100), hp=300))

#debugonce(ggplot2:::ggplot_build)

print(z1)
