library(tdmore)
library(nlmixr)
library(ggplot2)
library(testthat)

# Load the default tdmore
m1 <- getModel("default")
regimen <- data.frame(
  TIME=seq(0, 1)*24,
  AMT=5 #5mg
)

test_that("negative predictions should not mess up the error model", {
  m2 <- RxODE({
    TVKA = 3.7
    TVV1 = 61
    TVQ = 10
    TVCL = 3.7
    EPS_PROP = 0.23
    KA = TVKA * exp(EKA)
    CL = TVCL * (WT/70)^0.75 * exp(ECL)
    V1 = TVV1 * exp(EV1)
    V2 = V1
    Q = TVQ
    K12 = Q/V1
    K21 = Q/V2
    d/dt(depot) = -KA * depot
    d/dt(center) = KA * depot - CL/V1 * center - K12 * center +
      K21 * periph
    d/dt(periph) = K12 * center - K21 * periph
    CONC = center/V1
  })  %>% tdmore::tdmore(
    omega=c(EKA=0.4, ECL=0.0784, EV1=0.0361),
    res_var=list(errorModel(var="depot", prop=0.30))
  )
  tdmore:::predict.tdmore(m2, regimen=regimen, covariates=c(WT=70), newdata=seq(0, 48))
  expect_failure( #no warnings should be produced
    expect_warning(
      tdmore:::estimate(m2, regimen=regimen, covariates=c(WT=70), observed=data.frame(TIME=c(0.5, 1), depot=c(0.5, 0.02)))
    )
  )
})

test_that("Test all other ebe methods", {
  observed <- data.frame(TIME=10, CONC=0.04)
  fit <- tdmore:::estimate(m1, observed=observed,
                  regimen=regimen, covariates=c(WT=70))
  expect_snapshot_output( print(fit) )
  expect_snapshot_output( summary(fit) )
  expect_snapshot_value( confint(fit), style="serialize", tolerance=1E-4 )
  expect_equal( tdmore:::model.frame.tdmorefit(fit), observed )
  expect_snapshot_output( print(fitted(fit)) )
  expect_equal( logLik(fit), 3.98965892617313 )
  a <- logLik(fit, type="pop" )
  b <- logLik(fit, type="pred" )
  expect_equal( logLik(fit), a+b )
  expect_equal( a, 0.682675762497937 )
  expect_equal( b, 3.30698316367519 )
  expect_error( logLik(fit, type="babla" ) )
  expect_true( is.tdmorefit(fit) )
  expect_true( !is.tdmorefit(m1) )

  #How would a subject with higher bodyweight look? Can we re-estimate?
  fit2 <- tdmore:::estimate(fit, covariates=c(WT=75))
  expect_snapshot_output( summary(fit) )

  # test predict
  expect_equal( predict(fit), data.frame(TIME=10, CONC=0.03142424) , tolerance=1E-6)
  # test predict with fixed values
  par <- coef(fit)
  expect_equal( predict(fit, parameters=par[1]*3 ), data.frame(TIME=10, CONC=0.03550358), tolerance=1E-6)
  # test predict with fixed values and se.fit
  predict(fit, se.fit=T)
  predict(fit, se.fit=T, parameters=par[1]*3)
  predict(fit, se.fit=T, level=NA)
  foo <- predict(fit, se.fit=T, level=NA, parameters=par[1]*3)
  expect_equal(foo$ECL, rep(-0.5314563, 100), tolerance=1E-7 )
})

test_that("We can use multiple methods and pick the best one", {
  observed <- data.frame(TIME=10, CONC=0.04)
  fit <- tdmore:::estimate(m1, observed=observed, regimen=regimen, covariates=c(WT=70))

  cg <- tdmore:::estimate(fit, lower=NULL, upper=NULL, method="CG")
  expect_equal( coef(cg), coef(fit), tolerance=1e-6 )

  assign("deps", 1e-4, optextras::optsp) #fix for grfwd in optextras
  rcgmin <- tdmore:::estimate(fit, lower=NULL, upper=NULL, method="Rcgmin")
  nlminb <- tdmore:::estimate(fit, lower=NULL, upper=NULL, method="nlminb")

  #SANN returns all 0, because it does not work
  #And it takes a long time to boot...
  #sann <- tdmore:::estimate(fit, lower=NULL, upper=NULL, method="SANN")
  #expect_true( all(coef(sann) == 0) )

  multi <- tdmore:::estimate(fit, lower=NULL, upper=NULL, method=c("nlminb", "Rcgmin"))
  expect_equal( coef(multi), coef(rcgmin) )
  expect_true( logLik(rcgmin) > logLik(nlminb) )

  #allmeth <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", #"nlm", #nlm has errors with Grad
  #             "nlminb",
  #             "Rcgmin", "Rvmmin", "hjn")


  out <- capture.output(
    multistart <- tdmore:::estimate(fit, multistart=TRUE)
  )

  out <- capture.output(
    multistart <- tdmore:::estimate(fit, multistart=c(ECL=2, EV1=2))
  )
  out <- capture.output(
    multistart <- tdmore:::estimate(fit, multistart=c(ECL=2))
  )

  fitNoSe <- tdmore:::estimate(fit, se.fit=F)
  expect_true( all( diag(vcov(fitNoSe)) < 1e-8 ) )

  fixed <- tdmore:::estimate(fit, fix=c(EV1=1.2))
  expect_equal( coef(fixed)['EV1'] %>% unname , 1.2 )

  out <- capture.output(
    multiFix <- tdmore:::estimate(fit, fix=c(EV1=1.2), multistart=c(ECL=2))
  )
  expect_equal(coef(multiFix), coef(fixed), tolerance=1e-5)
  expect_equal(logLik(multiFix), logLik(fixed))

  expect_error({
    tdmore:::estimate(fit, fix=c(EV1=1, EV1=2))
  }) # EV1 appears twice!
  expect_error({
    tdmore:::estimate(fit, fix=c(EV1=1, ECL=2))
  }) #all parameters fixed...
})

describe("sampleMC is aware of IOV", {
  m2 <- getModel("example_iov")
  regimen <- data.frame(TIME=seq(0, by=8, length.out=12),
                        OCC=rep(1:4, each=3),
                        AMT=5)
  tdmore:::predict.tdmore(m2, regimen, newdata=c(0, 1, 2, 3)*24)

  ## Regimen in the tdmorefit object
  population <- tdmore:::estimate(m2, regimen=regimen)
  ipred <- predict(population, newdata=c(0, 1, 2, 3)*24, se.fit=TRUE, mc.maxpts=1, level=NA)
  expect_equal( unique(ipred$CL) %>% length , 4)
  ECL_IOV_Columns <- which( colnames(ipred) == "ECL_IOV")
  expect_equal( ipred$CL, 9.87 * exp(ipred$ECL + diag(as.matrix(ipred[,ECL_IOV_Columns]))))

  ## Extra occasion in the predict call
  ipred <- predict(population,
                   regimen=bind_rows(
                     regimen, tibble(TIME=100, AMT=10, OCC=5)
                     ), newdata=c(0, 24, 48, 72, 110)+5, se.fit=TRUE, mc.maxpts=1, level=NA)
  expect_equal( unique(ipred$CL) %>% length , 5)
  ECL_IOV_Columns <- which( colnames(ipred) == "ECL_IOV")
  expect_equal( length(ECL_IOV_Columns), 5)
  expect_equal( ipred$CL, 9.87 * exp(ipred$ECL + diag(as.matrix(ipred[,ECL_IOV_Columns]))))

  ## No occasions in the original tdmorefit
  expect_error(tdmore:::estimate(m2), "This model has IOV, a regimen has to be specified")
})


describe("We can generate uncertainty using MCMC", {
  observed <- data.frame(TIME=c(10, 12), CONC=c(0.04, 0.03))
  tdmorefit <- tdmore:::estimate(m1, observed=observed, regimen=regimen, covariates=c(WT=70))

  N <- 5000
  out1 <- sampleMC_norm(tdmorefit, mc.maxpts=N)$mc
  expect_equal(
    mean(out1$ECL),
    coef(tdmorefit)['ECL'],
    tolerance=0.05,
    ignore_attr=TRUE
  ) # should have the same mean

  set.seed(1234)
  out2 <- sampleMC_metrop(tdmorefit, mc.maxpts=N)$mc
  if(interactive()) {
    ggplot(out2 %>% tidyr::pivot_longer(cols=c(ECL, EV1)), aes(x=chain.sample, y=value, color=factor(chain), group=chain)) +
      geom_line() +
      facet_wrap(~name)
  }
  it("has a mean around the mode", {
    expect_equal(
      mean(out2$ECL),
      coef(tdmorefit)['ECL'],
      tolerance=0.2,
      ignore_attr=TRUE
    ) # should have the same mean, more or less
  })
  it("generates sufficient samples", {
    expect_gte(nrow(out2), N)
  })

  it("looks like fuzzy catterpillars", {
    p <- c(0.1, 0.5, 0.9)
    p_funs <- purrr::map(p, ~purrr::partial(stats::quantile, probs = .x, na.rm = TRUE)) %>%
      rlang::set_names(p)
    out2$slice <- cut(out2$chain.sample, breaks=seq(0, 10000, by=100) )
    summary <- out2 %>%
      dplyr::group_by(chain, slice) %>%
      dplyr::summarize_at(vars(ECL), p_funs) %>%
      tidyr::pivot_longer(cols = 3:5) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(value)
    q10 <- which(summary$name == "0.1")
    q50 <- which(summary$name == "0.5")
    q90 <- which(summary$name == "0.9")
    howManyBelow <- mean(q10 < min(q50) & q10 < min(q90 ) )
    expect_gte(howManyBelow, 0.9) #90% of 0.1 quantiles should be below all the others
    howManyAbove <- mean(q90 > max(q50) & q90 > max(q10) )
    expect_gte(howManyBelow, 0.9) #90% of 0.9 quantiles should be above all the others
  })

  ggplot(out2, aes(x=chain.sample, y=EV1, color=factor(chain))) + geom_step() +
    geom_hline(data=. %>% group_by(chain) %>% summarize(yintercept=mean(EV1)), aes(yintercept=yintercept))

  ggplot() + aes(x=EV1) + stat_density(data=out1, aes(color="rnorm"), fill=NA) + stat_density(data=out2, aes(color="MC"), fill=NA)
  ggplot() + aes(x=ECL) + stat_density(data=out1, aes(color="rnorm"), fill=NA) + stat_density(data=out2, aes(color="MC"), fill=NA)

  ggplot(mapping=aes(x=EV1, y=ECL))+
    stat_density2d(data=out1, aes(color="rnorm")) +
    stat_density2d(data=out2, aes(color="MC"))

  autoplot(tdmorefit, newdata=seq(0, 24, by=0.1))
  tdmorefit$varcov <- NULL
  autoplot(tdmorefit, newdata=seq(0, 24, by=0.1))
})
