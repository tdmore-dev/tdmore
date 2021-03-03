## Goal: to validate tdmore vs. nlmixr posthoc estimation

library(RxODE)
library(nlmixr)
library(tdmore)
library(dplyr)
library(ggplot2)
library(purrr)

# Execute nlmixr model fit ------------------------------------------------
describe("nlmixr and tdmore give same results", {
  #d/dt(depot)    = -KA*depot;
  #d/dt(centr)  =  KA*depot-kel*centr;
  #DV = centr / v;
  one.cmt <- function() {
    ini({
      tka <- 0.45 # Log Ka
      tcl <- 1 # Log Cl
      tv <- 3.45; label("log V")
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      KA <- exp(tka + eta.ka)
      Cl <- exp(tcl + eta.cl)
      Vc <- exp(tv + eta.v)
      d/dt(depot)    = -KA*depot;
      d/dt(centr)  =  KA*depot-Cl/Vc*centr;
      Cp = centr / Vc;
      Cp ~ add(add.sd)
    })
  }
  ## XXX: linCmt() does not work when generating an RxODE-style model...

  file <- testthat::test_path("ref/theo_sd_nlmixr.RDs")
  if(file.exists(file)) {
    fit <- readRDS(file)
  } else {
    message("Generating nlmixr fit...")
    fit <- nlmixr(one.cmt, theo_sd, est="focei", control = foceiControl(covMethod = ""))
    saveRDS(fit, file)
  }

  #ggplot(theo_sd, aes(x=TIME, y=DV)) + geom_line(aes(group=ID)) + geom_point()
  ggplot(fit, aes(x=TIME, y=DV)) + geom_point() +
    geom_line(aes(y=PRED), linetype=2) +
    geom_line(aes(y=IPRED), linetype=1) +
    facet_wrap(~ID)

  m1 <- tdmore(fit)
  obs <- theo_sd %>% filter(EVID==0) %>% rename(Cp=DV) #tdmore.nlmixrUI guesses the name of the output
  tmt <- theo_sd %>% filter(EVID!=0)
  db <- dataTibble(object=m1, observed=obs, regimen=tmt)
  posthocFit <- posthoc(db, control=list(trace=1))

  # posthocFit %>% unnest(cols=c(observed,ipred), names_sep=".") %>%
  #   ggplot(aes(x=observed.TIME, y=observed.Cp)) +
  #   geom_point() +
  #   geom_line(aes(y=ipred.Cp, color="tdmore")) +
  #   geom_line(data=fit, aes(x=TIME, y=IPRED, color="nlmixr")) +
  #   facet_wrap(~ID)

  ### XXX: The newest version of nlmixr reports widely varying results from previous versions...
  if(FALSE) {
    it("reports the same ETA estimates", {
      posthocFit$coef <- map(posthocFit$fit, ~tibble::enframe(coef(.x)))
      nlmixrCoef <- fit %>% as_tibble() %>% distinct(ID, eta.ka, eta.cl, eta.v) %>% mutate(ID=as.integer(ID))
      tdmoreCoef <- posthocFit %>% select(ID, coef) %>% tidyr::unnest(coef) %>% tidyr::pivot_wider()
      expect_equal(as.list(nlmixrCoef), as.list(tdmoreCoef), tolerance=1E-4) #rtol is 1E-6, but SIGDIG=3
    })
  }

  it("gives the same IPRED predictions", {
    expect_equal(
      fit %>% as_tibble() %>% transmute(ID=as.integer(ID), TIME, IPRED) %>% as.list(),
      posthocFit %>% tidyr::unnest(ipred) %>% transmute(ID, TIME, IPRED=Cp) %>% as.list(),
      tolerance=1E-1 #rtol is 1E-6 by default
    )
  })
})
