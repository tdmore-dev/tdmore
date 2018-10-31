library(nlmixr)
library(tdmore)
library(testthat)
library(RxODE)

context("Test the predict function (on both tdmore and tdmorefit objects)")

source(paste0(test_path(), ("/modelLibrary.R")))
m1 <- nlmixrUI(fluticasone_model) %>% tdmore()

regimen <- data.frame(
  TIME=seq(0, by=24, length.out=10),
  AMT=500 # 500ug standard dose
)

plot(m1, regimen, newdata = data.frame(TIME=seq(0, 250, by=1), CONC=NA))

## Let tdmore figure it out
## TDMore cannot figure it out, because we have no observations)
expect_error(predict(m1, regimen=regimen, newdata=NULL))  #Error: newdata is not a numeric or integer vector, indeed, newdata is a required argument

## You get all the variables
data <- predict(m1, regimen=regimen, newdata=seq(0, 30))
allVariables <- c("time", "Ka", "CL", "V1", "V2", "Q", "K12", "K21", "CONC", "center", "periph", "TIME")
expect_equal(colnames(data), allVariables)

## This is an error!
expect_error(predict(m1, regimen=regimen, newdata=data.frame(TIME=seq(0, 30)))) #Error: No output variable defined in newdata

pred <- estimate(m1, regimen=regimen)
ipred <- estimate(m1, observed=data.frame(TIME=75, CONC=22), regimen=regimen)

## TDMore guesses what you want
## This behaviour is very badly defined....
## Actually, we should simply show the prediction for the 'observed' data
## If observed data is empty, we should give back an empty data.frame
data <- predict(pred, regimen=regimen, newdata=NULL, se.fit=FALSE)
expect_true(is.data.frame(data) && nrow(data)==0)

data <- predict(ipred, regimen=regimen, newdata=NULL, se.fit=FALSE)
expect_equal(round(data, 2), data.frame(TIME=75, CONC=21.92))

## Great, we get all the variables
data <- predict(pred,
        regimen = regimen,
        newdata = seq(0, 30),
        se.fit = FALSE)
expect_equal(colnames(data), allVariables)

## We should get all of the parameters
data <- predict(pred,
                regimen = regimen,
                newdata = seq(0, 30),
                se.fit = TRUE)
# Check only the first five columns
expect_equal(colnames(data)[1:5], c("TIME", "Ka", "Ka.median", "Ka.lower", "Ka.upper"))

## Ok, also an error, same as for predict.tdmore
expect_error(predict(pred, regimen=regimen,
        newdata=data.frame(TIME=seq(0, 30)),
        se.fit=TRUE)) # Error: No output variable defined in newdata

## Now we get TIME and CONC
data <- predict(pred, regimen=regimen,
        newdata=data.frame(TIME=seq(0, 30), CONC=NA),
        se.fit=TRUE)
expect_equal(colnames(data), c("TIME", "CONC", "CONC.median", "CONC.lower", "CONC.upper"))
