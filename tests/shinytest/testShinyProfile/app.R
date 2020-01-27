library(tdmore)
library(RxODE)
library(nlmixr)
library(ggplot2)

### ensure input/output/export all have a consistent order
### Sorting order is different between en_US and C locales!
### See ?Comparison
### Fix is to force the locale to C for collation
Sys.setlocale(category="LC_COLLATE", "C")

tdmore <- getModel("meropenem")
regimen <- data.frame(
  TIME=c(0, 8, 16),            # Every 8 hour and for 1 day, an injection is given
  AMT=c(1000, 1000, 1000),     # 1g is administered
  RATE=c(1000, 1000, 1000)/0.5 # 30-minute infusion (rate=dose/infusion time)
)
data <- predict(
  object = tdmore,
  newdata = data.frame(TIME = seq(0, 24, by = 0.5), CONC = NA),
  regimen = regimen,
  se = TRUE
)

pred <- tdmore %>% estimate(regimen = regimen)
observed <- data.frame(TIME=c(9, 16), CONC=c(30, 8))
ipred <- tdmore %>% estimate(observed = observed, regimen = regimen)

shinyProfileApp(ipred)
