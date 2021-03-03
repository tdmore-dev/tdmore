library(tdmore)
library(RxODE)
library(nlmixr)
library(ggplot2)

### ensure input/output/export all have a consistent order
### Sorting order is different between en_US and C locales!
### See ?Comparison
### Fix is to force the locale to C for collation
Sys.setlocale(category="LC_COLLATE", "C")

model <- tdmore:::getModel("example")
regimen <- data.frame(
  TIME=c(0, 8, 16),            # Every 8 hour and for 1 day, an injection is given
  AMT=c(1000, 1000, 1000),     # 1g is administered
  RATE=c(1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
)

observed <- data.frame(TIME=c(9, 16), CONC=c(30, 8))
ipred <- model %>% tdmore:::estimate(observed = observed, regimen = regimen)
ipred$res <- round(ipred$res, digits=6) #round so inter-machine differences do not pop up...
ipred$varcov <- round(ipred$varcov, digits=6)

app <- tdmore:::shinyProfileApp(ipred)

app
