modelCode <- "
CL = 23.6 * exp(ETA1*0.42);
Vc = 1070 * exp(ETA2*1.11);
ka=4.48;
CONC = centr / Vc * 1000;
d/dt(abs) = -ka*abs;
d/dt(centr) = ka*abs - CL/Vc*centr;
"
rxodeModel <- RxODE::RxODE(modelCode)
model <- tdmore(rxodeModel, res_var=list(errorModel(var="CONC", add=3.7)))
