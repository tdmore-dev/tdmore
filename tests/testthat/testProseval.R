
m1 <- getModel("default")
regimen <- data.frame(TIME=0, AMT=1000, DURATION=0.5)
cov <- c(WT=70)

describe("posthoc", {
  it("returns estimates for multiple subjects", {
    observed1 <- data.frame(ID=1, TIME=3, CONC=8)
    fit1 <- estimate(m1, observed=observed1, regimen=regimen, covariates=cov)

    observed2 <- data.frame(ID=2, TIME=3, CONC=12)
    fit2 <- estimate(m1, observed=observed2, regimen=regimen, covariates=cov)

    db <- dataTibble(object=m1, observed=rbind(observed1, observed2))
    expect_error({posthoc(db, regimen=regimen)}, "WT is missing")
    fits <- posthoc(db, regimen=regimen, covariates=cov)

    expect_equal( #columns from db are repeated
      db %>% as.list(),
      fits %>% dplyr::select(!!colnames(db)) %>% as.list(),
      ignore_attr=TRUE
    )
    expect_equal(
      setdiff(names(fits), colnames(db)),
      c("fit", "elapsed", "ipred") #other columns: fit, elapsed time, ipred
    )

    expect_equal(fits$ID, c(1,2))
    expect_equal( coef(fits$fit[[1]]),
                  coef(fit1)
    )
    expect_equal( coef(fits$fit[[2]]),
                  coef(fit2)
    )
  })
})

describe("proseval", {
  it("splits a concentration profile into multiple sub-parts", {
    observed <- data.frame(TIME=c(3, 6), CONC=c(8, 4))
    fit0 <- estimate(m1, observed=observed[0,,drop=F], regimen=regimen, covariates=cov)
    fit1 <- estimate(m1, observed=observed[1,,drop=F], regimen=regimen, covariates=cov)
    fit2 <- estimate(m1, observed=observed[1:2,,drop=F], regimen=regimen, covariates=cov)
    reference <- list(fit0, fit1, fit2)

    db <- dataTibble(object=m1, observed=observed)
    fits <- proseval(db, regimen=regimen, covariates=cov)

    for(i in seq_along(reference)) {
      expect_equal( coef(reference[[i]]),
                    coef(fits$fit[[i]])
      )
    }
  })
  it("uses the par argument, also if it is a list", {
    observed <- data.frame(TIME=c(3, 6), CONC=c(8, 4))
    db <- dataTibble(object=m1, observed=observed)
    fitsOutput <- capture.output(
      {proseval(db, regimen=regimen, covariates=cov, control=list(trace=2))}
    )

    fitsOutput2 <- capture.output(
      {proseval(db, par=c(ECL=2, EV1=2), regimen=regimen, covariates=cov, control=list(trace=2))}
    )

    db$par <- NA
    db$par[1] <- list( lapply(list(fit1, fit2), coef) )
    fitsOutput3 <- capture.output(
      {proseval(db, regimen=regimen, covariates=cov, control=list(trace=2))}
    )

    for(i in seq_along(reference)) {
      expect_equal( coef(reference[[i]]),
                    coef(fits$fit[[i]])
      )
    }
  })

  it("should *not* split the covariates as well", {
    observed <- data.frame(TIME=c(3, 6), CONC=c(8, 4))
    covariates <- tibble(TIME=c(0, 2, 5, 9), WT=c(70, 72, 74, 76))
    fit0 <- estimate(m1, observed=observed[0,,drop=F], regimen=regimen, covariates=covariates)
    fit1 <- estimate(m1, observed=observed[1,,drop=F], regimen=regimen, covariates=covariates)
    fit2 <- estimate(m1, observed=observed[1:2,,drop=F], regimen=regimen, covariates=covariates)
    reference <- list(fit0, fit1, fit2)

    db <- dataTibble(object=m1, observed=observed)
    fits <- proseval(db, regimen=regimen, covariates=covariates)

    for(i in seq_along(reference)) {
      expect_equal( coef(reference[[i]]),
                    coef(fits$fit[[i]])
      )
    }
  })
})

describe("doseSimulation", {
  it("Simulates a dose adaptation", {
    set.seed(1234)
    observed <- data.frame(TIME=c(3, 6), CONC=c(8, 4))
    regimen <- data.frame(TIME=seq(0, 24, by=8), AMT=1000, DURATION=0.5)
    db <- dataTibble(object=m1, observed, regimen, covariates=cov)
    fits <- posthoc(db) %>% dplyr::select(-elapsed, -ipred)
    sim <- doseSimulation(fits, control=list(trace=1), optimize=function( fit, regimen, truth ){
      rows <- which(regimen$TIME > max(c(-Inf, fit$observed$TIME)))
      if(length(rows)==0) return(regimen) #nothing to adapt
      rec <- findDose(fit, regimen, doseRows=rows, target=data.frame(TIME=24+8, CONC=16))
      list( regimen=rec$regimen, extra=list(recommendation=list(rec)) )
    })

    z <- autoplot(sim$iterationFit[[1]], newdata=seq(0, 32, by=0.1)) +
      autolayer(sim$recommendation[[1]]) +
      autolayer(sim$fit[[1]], se.fit=F, regimen=sim$recommendation[[1]]$regimen, linetype=2, color="purple") +
      coord_cartesian(ylim=c(0, 50))
    if(interactive()) expect_doppelganger("Recommendation hits target_1", z)

    z <- autoplot(sim$iterationFit[[2]], newdata=seq(0, 32, by=0.1)) +
      autolayer(sim$recommendation[[2]]) +
      autolayer(sim$fit[[1]], se.fit=F, regimen=sim$recommendation[[2]]$regimen, linetype=2, color="purple") +
      coord_cartesian(ylim=c(0, 50))
    if(interactive()) expect_doppelganger("Recommendation hits target_2", z)

    z <- autoplot(sim$iterationFit[[3]], newdata=seq(0, 32, by=0.1)) +
      autolayer(sim$recommendation[[3]]) +
      autolayer(sim$fit[[1]], se.fit=F, regimen=sim$recommendation[[3]]$regimen, linetype=2, color="purple") +
      coord_cartesian(ylim=c(0, 50))
    if(interactive()) expect_doppelganger("Recommendation hits target_3", z)

  })
})
