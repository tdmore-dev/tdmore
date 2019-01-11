library(tdmore)


## We can now repeat the exercise for a different patient
covariates <- c(BW=70, CRCL=56, AGE=73)
regimen <- bind_rows(
  data.frame(TIME=0, II=12, ADDL=1, AMT=750),
  data.frame(TIME=24, II=12, ADDL=11, AMT=1000),
  data.frame(TIME=7*24, II=12, ADDL=5, AMT=500)
  #TDMore will always try to adapt the last row in the regimen
)
observed <- bind_rows(
  data.frame(TIME=24, CONC=12.5),
  data.frame(TIME=3*24, CONC=17.5),
  data.frame(TIME=6*24, CONC=20.2),
  data.frame(TIME=7*24, CONC=17.1),
  data.frame(TIME=9*24, CONC=15.2)
)
plot(model, regimen=regimen, covariates=covariates, newdata=seq(0, 10*24, by=0.1)) +
  geom_point(aes(y=CONC), data=observed) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  scale_x_continuous(breaks=seq(0, 10*24, by=12))

ipred <- estimate(model, regimen=regimen, covariates=covariates, observed=observed[1,])
plot(ipred) +
  geom_hline(yintercept=c(TARGET, SAFETY)) +
  geom_point(data=observed, aes(y=CONC)) +
  scale_x_continuous(breaks=seq(0, 10*24, by=12))
