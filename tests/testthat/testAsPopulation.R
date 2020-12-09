library(dplyr)

test_that("as.sample", {
  m1 <- getModel()
  set.seed(1234)
  db <- m1 %>%
    as.population(covariates=c(WT=70)) %>%
    as.sample(N=10)

  ## Plot result with 5mg dosing
  predictions <- db %>% mutate(
    ipred = purrr::map(fit, predict, regimen=data.frame(TIME=0, AMT=5), newdata=0:24)
  )
  z1 <- predictions %>% tidyr::unnest(cols=ipred) %>%
    ggplot(aes(x=TIME, y=CONC)) +
    geom_line(aes(group=ID))
  expect_doppelganger("population-sample", z1)

  ## Plot result with optimized dosing
  optimResult <- doseSimulation(db, regimen=data.frame(TIME=0, AMT=5),
    optimize=function(fit, regimen, truth) {
      rec <- findDose(fit, target = data.frame(TIME=24, CONC=0.05))
      list(nextTime=if(nrow(fit$observed)==0) 12 else NA,
           regimen=rec$regimen)
  })
  predictions <- optimResult %>% ungroup() %>% mutate(
    ipred = purrr::map2(fit, next_regimen, ~predict(.x, regimen=.y, covariates=c(WT=70), newdata=0:24))
  )
  z1 <- predictions %>% tidyr::unnest(cols=ipred) %>%
    ggplot(aes(x=TIME, y=CONC, color=factor(OBS))) +
    geom_line(aes(group=interaction(ID,OBS)))
  expect_doppelganger("doseSimulation-sample", z1)
})
