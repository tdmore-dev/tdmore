context("Proseval and other multi-ID methods work correctly")

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

    expect_equivalent( #columns from db are repeated
      db %>% as.list(),
      fits %>% select(!!colnames(db)) %>% as.list()
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
})

describe("doseSimulation", {
  it("Simulates a dose adaptation", {
    observed <- data.frame(TIME=c(3, 6), CONC=c(8, 4))
    regimen <- data.frame(TIME=seq(0, 24, by=8), AMT=1000, DURATION=0.5)
    db <- dataTibble(object=m1, observed, regimen, covariates=cov)
    fits <- posthoc(db) %>% select(-elapsed, -ipred)
    sim <- doseSimulation(fits, control=list(trace=1), optimize=function( fit, regimen, truth ){
      rows <- which(regimen$TIME > max(c(-Inf, fit$observed$TIME)))
      if(length(rows)==0) return(regimen) #nothing to adapt
      rec <- findDose(fit, regimen, doseRows=rows, target=data.frame(TIME=24+8, CONC=16))
      list( regimen=rec$regimen, extra=list(recommendation=rec) )
    })

    z <- autoplot(sim$fit[[1]], newdata=seq(0, 32, by=0.1)) +
      autolayer(sim$recommendation[[1]])
    expect_doppelganger("Recommendation hits target", z)

  })
})
