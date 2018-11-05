library(RxODE)

mod1 <-RxODE({
  C2 = centr/V2;
  C3 = peri/V3;
  d/dt(depot) =-KA*depot;
  d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
  d/dt(peri)  =                    Q*C2 - Q*C3;
  d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
});

theta <-
  c(KA=2.94E-01, CL=1.86E+01, V2=4.02E+01, # central
    Q=1.05E+01,  V3=2.97E+02,              # peripheral
    Kin=1, Kout=1, EC50=200)               # effects

ev <- eventTable(amount.units='mg', time.units='hours')
ev$add.dosing(dose=10000, nbr.doses=280, dosing.interval=24, start.time=0)

#times <- seq(0, 40*7*24,by=1) # TRY WITH THIS TIMES VECTOR: OK
times <- c(0, seq(3990, 4010, by = 1), 40*7*24) # TRY WITH THIS ONE: NOK -> 5000 steps taken before reaching tout

ev$add.sampling(times)

knitr::kable(head(ev$get.dosing()))
knitr::kable(head(ev$get.sampling()))

inits <- c(eff=1)

x <- solve(mod1, theta, ev, inits);
x <- as.data.frame(x)
ggplot(x,aes(time,C2)) + geom_line() + ylab("Central Concentration") + xlab("Time")
