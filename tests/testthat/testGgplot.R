library(tdmore)
library(nlmixr)
library(ggplot2)
library(testthat)

context("Test ggplot functions")

m1 <- tdmore(pheno_nlmixr)
regimen <- data.frame(
  TIME=seq(0, 1)*24,
  AMT=5 #5mg
)
observed <- data.frame(TIME=c(2, 5), CONC=c(0.040, 0.02))
pred <- m1 %>% estimate(regimen = regimen, covariates=c(WT=70))
ipred <- m1 %>% estimate(observed = observed, regimen = regimen, covariates = c(WT=70))

test_that("stat_predict() works like a cross-breed of predict() and stat_function()", {
  predict(ipred, newdata=0:10)

  # some examples from stat_function
  ggplot() + stat_function(fun=sin) #empty, no data
  ggplot() + stat_function(fun=sin, xlim=c(0, 10)) #empty too, no axis was defined
  ggplot() + stat_function(fun=sin, data=data.frame(x=0:10), mapping=aes(x)) #data!
  ggplot(data.frame(x=0:10)) + stat_function(fun=sin, mapping=aes(x)) #data!
  #expect_warning(  #weirdly enough, this says "code did not produce any warnings"
  #  ggplot(data.frame(x=0:10)) + stat_function(fun=sin) #no data on x axis :(
  #)
  # Or with xlim:
  ggplot() +
    geom_histogram(data=data.frame(x=runif(100)), aes(x)) +
    stat_function(fun=sin, xlim=c(0,1), color="red") #no line
  ggplot(data=data.frame(x=runif(100))) +
    geom_histogram(aes(x)) +
    stat_function(fun=sin, xlim=c(0,1), color="red") #line! because common axis was defined?

  #this is the goal
  ggplot() + stat_function(fun=function(x, tdmorefit){
      result <- predict(tdmorefit, newdata=x)
      result$CONC
    }, args=list(tdmorefit=ipred), data=data.frame(x=0:50), mapping=aes(x))

  # using stat_predict
  ggplot() +
    stat_predict(tdmorefit=ipred, data=data.frame(TIME=0:50), mapping=aes(x=TIME))
  # Moving stuff to top ggplot
  ggplot(ipred) +
    stat_predict(data=data.frame(TIME=0:50), mapping=aes(x=TIME))

  # Move data to top
  ggplot(ipred, aes(x=TIME)) +
    geom_point(aes(y=CONC, color="data")) +
    stat_predict(mapping=aes(x=TIME, color="prediction"), geom="point")

  ggplot(ipred, aes(x=TIME)) +
    geom_point(aes(y=CONC, color="data")) +
    stat_predict(mapping=aes(x=TIME, color="prediction"), geom="point", xlim=c(NA,NA))

  ggplot(ipred) +
    stat_predict(mapping=aes(x=TIME), xlim=c(0, 12), geom="point")

  ggplot(ipred, aes(x=TIME, y=CONC)) +
    geom_point() +
    stat_predict(mapping=aes(x=TIME), xlim=c(0, 12))

  ggplot(ipred, aes(x=TIME, y=CONC)) +
    geom_point() +
    stat_predict(mapping=aes(x=TIME), xlim=c(0, 12))

  # This means the data can be moved to top as well
  ggplot(pred, newdata=data.frame(TIME=0:50)) +
    stat_predict(mapping=aes(x=TIME))

  ## X should be available for stat_predict
  z1 <- ggplot() + stat_predict(tdmorefit=pred, data=data.frame(TIME=0:50))
  expect_warning( layer_data(z1) ) #plot needs to build to get the error
})

foo <- predict(object=pred, newdata=seq(0, 50))
ggplot(foo) + stat_identity(aes(x=TIME, y=CONC))

z1 <- ggplot(pred, aes(x=TIME,y=CONC)) +
  stat_predict(geom="line",
               xlim=c(0, 10))
print(z1)


########################################

# library(tdmore)
# library(ggplot2)
# m1 <- tdmore(meropenem_model)
# regimen <- data.frame(AMT=2000, TIME=seq(0, 100, by=8))
#
# population <- estimate(m1, regimen=regimen)
# observed <- data.frame(TIME=c(4, 7.8), CONC=c(10, 2.4))
# ipred <- estimate(m1, regimen=regimen, observed=observed)
#
# ggplot(ipred, mapping=aes(x=TIME,y=CONC)) +
#   geom_point() +
#   stat_predict(xlim=c(NA,NA)) +
#   stat_predict(xlim=c(NA,NA), tdmorefit=population)
#
# ggplot(ipred, mapping=aes(x=TIME,y=CONC)) +
#   geom_point() +
#   geom_fit(xlim=c(NA,NA))
