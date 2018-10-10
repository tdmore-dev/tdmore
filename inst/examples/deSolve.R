# # Example for using deSolve with tdmore
# rm(list=ls(all=TRUE))
# func <- function(t, y, parms) {
#   ETA1 <- as.numeric(parms['ETA1'])
#   ETA2 <- as.numeric(parms['ETA2'])
#   CL = 23.6 * exp(ETA1*0.42)
#   Vc = 1070 * exp(ETA2*1.11);
#   ka=4.48;
#   abs = as.numeric(y[1])
#   centr = as.numeric(y[2])
#   CONC = centr / Vc * 1000
#   dabs = -ka*abs
#   dcentr = ka*abs - CL/Vc*centr
#   return(list(c(dabs,dcentr), c(CONC=CONC)))
# }
# model <- tdmore_deSolve(parameters=c("ETA1", "ETA2"),
#                         add=3.7,
#                         func=func,
#                         y=c(abs=0, centr=0),
#                         ynames=FALSE)
# regimen <- data.frame(
#   TIME=c(0, 24),
#   AMT=c(15, 15),
#   II=c(0, 24),
#   ADDL=c(0, 10)
# )
# pred <- predict(model, newdata=0:24, regimen=regimen)
# ggplot(pred, aes(x=TIME, y=CONC)) + geom_line()
#
# observed <- data.frame(TIME=c(2.4, 23), CONC=c(10, 5))
# ipredfit <- model %>%
#   estimate(observed, regimen)
# z <- plot(ipredfit)
# print(z)
#
# D <- findDose(ipredfit, regimen=regimen, target=data.frame(TIME=48, CONC=13.5))
