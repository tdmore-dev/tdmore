library(tdmore)
library(RxODE)
library(testthat)

context("Test that the algebraic library gives the same results as RxODE")

solveODE <- function(model, times, TIME, AMT, II=NULL, RATE=NULL, N=200, ...) {
  m1 <- RxODE(model)

  if(!is.null(II)) {
    ev <- eventTable() %>% add.dosing(dose=AMT, start.time=TIME, nbr.doses=N+1, dosing.interval=II, rate=RATE) %>%
      add.sampling(times + II*N)
  } else {
    ev <- eventTable() %>% add.dosing(dose=AMT, start.time=TIME, rate=RATE) %>%
      add.sampling(times)
  }

  params <- list(...)
  result <- rxSolve(m1, params = unlist(params), events = ev, atol=1E-12, rtol=1E-12, maxsteps=N*10*1000) %>% as.data.frame
  result$CONC
}

times <- seq(0, 24, by=0.1)
tD <- 0
AMT <- 50


# 1cpt models -------------------------------------------------------------
test_that("1cpt models", {
  model <- "
  d/dt(A1) = -K*A1;
  CONC = A1 / V;
  "
expect_equal(
  pk1cptiv_(times, tD, AMT, K=0.3, V=10),
  solveODE(model, times, tD, AMT, K=0.3, V=10)
)

expect_equal(
  pk1cptiv_(times, tD, AMT, K=0.03, V=20, II=24, SS=1),
  solveODE(model, times, tD, AMT, K=0.03, V=20, II=24)
)


expect_equal(
  pk1cptinfusion_(times, tD, AMT, K=0.3, V=10, RATE=25),
  solveODE(model, times, tD, AMT, K=0.3, V=10, RATE=25)
)

expect_equal(
  pk1cptinfusion_(times, tD, AMT, K=0.03, V=20, II=24, RATE=25, SS=1),
  solveODE(model, times, tD, AMT, K=0.03, V=20, II=24, RATE=25)
)
model <- "
d/dt(A0) = -KA*A0;
d/dt(A1) = KA*A0 -K*A1;
CONC = A1 / V;
"
expect_equal(
  pk1cptoral_(times, tD, AMT, KA=0.3, K=0.03, V=20),
  solveODE(model, times, tD, AMT, KA=0.3, K=0.03, V=20)
)

expect_equal(
  pk1cptoral_(times, tD, AMT, KA=0.3, K=0.03, V=20, II=24, SS=1),
  solveODE(model, times, tD, AMT, KA=0.3, K=0.03, V=20, II=24)
)
})


expect_warning(
pk1cptinfusion_(times, tD, AMT, K=0.03, V=20, II=24, RATE=2, SS=1),
"Infusion time larger than interdose interval"
)


# 2cpt models -------------------------------------------------------------
test_that("2cpt models", {
  model <- "
  d/dt(A1) = -K*A1 -K12*A1 + K21*A2;
  d/dt(A2) = K12*A1 - K21*A2;
  CONC = A1 / V;
  "
  expect_equal(
    pk2cptiv_(times, tD, AMT, K=0.3, V=10, K12=0.13, K21=0.08),
    solveODE(model, times, tD, AMT, K=0.3, V=10, K12=0.13, K21=0.08)
  )

  expect_equal(
    pk2cptiv_(times, tD, AMT, K=0.03, V=20, K12=0.13, K21=0.08, II=24, SS=1),
    solveODE(model, times, tD, AMT, K=0.03, V=20, K12=0.13, K21=0.08, II=24)
  )


  expect_equal(
    pk2cptinfusion_(times, tD, AMT, K=0.3, V=10, K12=0.13, K21=0.08, RATE=25),
    solveODE(model, times, tD, AMT, K=0.3, V=10, K12=0.13, K21=0.08, RATE=25)
  )

  expect_equal(
    pk2cptinfusion_(times, tD, AMT, K=0.03, V=20, K12=0.13, K21=0.08, II=24, RATE=25, SS=1),
    solveODE(model, times, tD, AMT, K=0.03, V=20, K12=0.13, K21=0.08, II=24, RATE=25)
  )
  model <- "
  d/dt(A0) = -KA*A0;
  d/dt(A1) = KA*A0 -K*A1 -K12*A1 + K21*A2;
  d/dt(A2) = K12*A1 - K21*A2;

  CONC = A1 / V;
  "
  expect_equal(
    pk2cptoral_(times, tD, AMT, KA=0.3, K=0.03, V=20, K12=0.13, K21=0.08),
    solveODE(model, times, tD, AMT, KA=0.3, K=0.03, V=20, K12=0.13, K21=0.08)
  )

  expect_equal(
    pk2cptoral_(times, tD, AMT, KA=0.3, K=0.03, V=20, K12=0.13, K21=0.08, II=24, SS=1),
    solveODE(model, times, tD, AMT, KA=0.3, K=0.03, V=20, K12=0.13, K21=0.08, II=24)
  )
})


# 3cpt --------------------------------------------------------------------
test_that("3cpt models", {
  model <- "
  d/dt(A1) = -K*A1 -K12*A1 + K21*A2 -K13*A1 + K31*A3;
  d/dt(A2) = K12*A1 - K21*A2;
  d/dt(A3) = K13*A1 - K31*A3;
  CONC = A1 / V;
  "
  expect_equal(
    pk3cptiv_(times, tD, AMT, K=0.3, V=10, K12=0.13, K21=0.08, K13=0.2, K31=0.09),
    solveODE(model, times, tD, AMT, K=0.3, V=10, K12=0.13, K21=0.08, K13=0.2, K31=0.09)
  )

  expect_equal(
    pk3cptiv_(times, tD, AMT, K=0.03, V=20, K12=0.13, K21=0.08, K13=0.2, K31=0.09, II=24, SS=1),
    solveODE(model, times, tD, AMT, K=0.03, V=20, K12=0.13, K21=0.08, K13=0.2, K31=0.09, II=24)
  )

  expect_equal(
    pk3cptinfusion_(times, tD, AMT, K=0.3, V=10, K12=0.13, K21=0.08, K13=0.2, K31=0.09, RATE=25),
    solveODE(model, times, tD, AMT, K=0.3, V=10, K12=0.13, K21=0.08, K13=0.2, K31=0.09, RATE=25)
  )

  expect_equal(
    pk3cptinfusion_(times, tD, AMT, K=0.03, V=20, K12=0.13, K21=0.08, K13=0.2, K31=0.09, II=24, RATE=25, SS=1),
    solveODE(model, times, tD, AMT, K=0.03, V=20, K12=0.13, K21=0.08, K13=0.2, K31=0.09, II=24, RATE=25)
  )
  model <- "
  d/dt(A0) = -KA*A0;
  d/dt(A1) = KA*A0 -K*A1 -K12*A1 + K21*A2 -K13*A1 + K31*A3;
  d/dt(A2) = K12*A1 - K21*A2;
  d/dt(A3) = K13*A1 - K31*A3;

  CONC = A1 / V;
  "
  expect_equal(
    pk3cptoral_(times, tD, AMT, KA=0.3, K=0.03, V=20, K12=0.13, K21=0.08, K13=0.2, K31=0.09),
    solveODE(model, times, tD, AMT, KA=0.3, K=0.03, V=20, K12=0.13, K21=0.08, K13=0.2, K31=0.09)
  )

  expect_equal(
    pk3cptoral_(times, tD, AMT, KA=3, K=0.03, V=20, K12=0.13, K21=0.08, K13=0.2, K31=0.9, II=24, SS=1),
    solveODE(model, times, tD, AMT, KA=3, K=0.03, V=20, K12=0.13, K21=0.08, K13=0.2, K31=0.9, II=24)
  )
})







# Test manual selection of model ---------------------------------------------------
m1 <- algebraic(function(t, TIME, AMT){
  KA=0.3
  K=0.03
  V=20
  K12=0.13
  K21=0.08
  K13=0.2
  K31=0.09
  II=24
  pk2cptiv()
})
model_predict(m1, times=seq(0, 24, by=0.1), regimen=data.frame(TIME=0, AMT=15))

m1 <- algebraic(function(t, TIME, AMT){
  KA=0.3
  K=0.03
  V=20
  K12=0.13
  K21=0.08
  K13=0.2
  K31=0.09
  II=24
  SS=1
  pk2cptoral()
})
model_predict(m1, times=seq(0, 24, by=0.1), regimen=data.frame(TIME=0, AMT=15))

expect_error({
  m1 <- algebraic(function(t, TIME, AMT){
    K=0.03
    V=20
    K12=0.13
    K21=0.08
    pk2cptoral()
  })
  model_predict(m1, times=seq(0, 24, by=0.1), regimen=data.frame(TIME=0, AMT=15))
}, ".*KA not found.*")

# Test automatic detection, nlmixr style ----------------------------------
m1 <- algebraic(function(t, TIME, AMT){
  KA=0.3
  KE=0.03
  V=20
  K12=0.13
  K21=0.08
  K13=0.2
  K31=0.09
  linCmt()
})
model_predict(m1, times=seq(0, 24, by=0.1), regimen=data.frame(TIME=0, AMT=15))

m1 <- algebraic(function(t, TIME, AMT){
  KA=0.3
  CL=0.03
  V=20
  CLD=0.13
  VT=0.08
  CLD2=0.2
  VT2=0.09
  linCmt()
})
model_predict(m1, times=seq(0, 24, by=0.1), regimen=data.frame(TIME=0, AMT=15))


m1 <- algebraic(function(t, TIME, AMT){
  KA=0.3
  CL=0.03
  linCmt()
})
expect_error(
model_predict(m1, times=seq(0, 24, by=0.1), regimen=data.frame(TIME=0, AMT=15)),
"No matching model.*")


# Test automatic detection, Monolix-style ---------------------------------
m1 <- algebraic(function(t, TIME, AMT){
  ka=0.3
  k=0.03
  V=20
  k12=0.13
  k21=0.08
  K13=0.2
  K31=0.09
  II=24
  pkmodel(ka, k, V, k12, k21)
})
model_predict(m1, times=seq(0, 24, by=0.1), regimen=data.frame(TIME=0, AMT=15))

m1 <- algebraic(function(t, TIME, AMT){
  V=20
  k12=0.13
  k21=0.08
  K13=0.2
  K31=0.09
  II=24
  pkmodel(ka=0.3, k=0.03, V, k12, k21)
})
model_predict(m1, times=seq(0, 24, by=0.1), regimen=data.frame(TIME=0, AMT=15))



m1 <- algebraic(function(t, TIME, AMT){
  V=-20
  k12=0.13
  k21=0.08
  K13=0.2
  K31=0.09
  II=24
  pkmodel(ka=0.3, k=0.03, V=-1*V+32, k12, k21)
})
model_predict(m1, times=seq(0, 24, by=0.1), regimen=data.frame(TIME=0, AMT=15))
