library(testthat)

context('Test model.frame functions')

m1 <- tdmore(meropenem_model)
test_that("model.frame correctly returns values", {
  expect_error( model.frame(m1) )
  expect_equal(
    model.frame(m1, data=0:10),
    data.frame(TIME=0:10, CONC=as.numeric(NA))
  )

  expect_equal(
    model.frame(m1, data=data.frame(TIME=0:10, CONC=13)),
    data.frame(TIME=0:10, CONC=13)
  )
})

test_that("model.frame correctly returns SE values", {
  expect_error( model.frame(m1, se=TRUE) )
  expect_equal(
    model.frame(m1, data=0:10, se=TRUE),
    data.frame(TIME=0:10, CONC=as.numeric(NA), CONC.lower=as.numeric(NA), CONC.upper=as.numeric(NA))
  )

  expect_equal(
    model.frame(m1, data=data.frame(TIME=10, CONC=20), se=TRUE),
    data.frame(TIME=10, CONC=20, CONC.lower=34.5429327652872, CONC.upper=5.4570672347128)
  )
})
