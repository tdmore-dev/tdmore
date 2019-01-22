library(testthat)
library(tdmore)
library(deSolve)

context("Test deSolve model")

func <- function(t, y, parms) {
  ETA1 <- as.numeric(parms['ETA1'])
  ETA2 <- as.numeric(parms['ETA2'])
  CL = 23.6 * exp(ETA1*0.42)
  Vc = 1070 * exp(ETA2*1.11);
  ka=4.48;
  abs = as.numeric(y[1])
  centr = as.numeric(y[2])
  CONC = centr / Vc * 1000
  dabs = -ka*abs
  dcentr = ka*abs - CL/Vc*centr
  return(list(c(dabs,dcentr), c(CONC=CONC)))
}
model <- tdmore_deSolve(parameters=c("ETA1", "ETA2"),
                        add=3.7,
                        func=func,
                        y=c(abs=0, centr=0),
                        ynames=FALSE)
regimen <- data.frame(
  TIME=c(0, 24),
  AMT=c(15, 15),
  II=c(0, 24),
  ADDL=c(0, 10)
)

pred <- predict(model, newdata=seq(0, 200, by=0.1), regimen=regimen)
expect_known_value(pred, "deSolvePrediction")
#ggplot(pred, aes(x=TIME, y=CONC)) + geom_line()
