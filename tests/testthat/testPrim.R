library(subgroup.discovery)
context("Patient Rule Induction Method")


testthat::test_that("Test the data preparation phase", {

  data(credit)
  noncompliant <- any(sapply(credit, function(x){!(is.numeric(x) | is.factor(x))}))
  credit <- prim.data.prepare(credit)
  compliant <- any(sapply(credit, function(x){!(is.numeric(x) | is.factor(x))}))

  expect_true(noncompliant)
  expect_false(compliant)
})

testthat::test_that("Test functionality on pima data set", {

  data(pima)
  pima.sample <- 1:(0.75*nrow(pima))
  pima <- prim.data.prepare(pima)
  pima.model <- prim(class ~ ., data = pima[pima.sample,], peeling.quantile = 0.4, min.support = 0.4)
  pima.predict <- predict(pima.model, pima[-pima.sample,])

  expect_is(pima.model, "prim.peel")
  expect_is(pima.predict, "prim.predict")

  pima.model.idx <- prim.box.index(pima.model, pima)
  pima.predict.idx <- prim.box.index(pima.predict, pima)

  expect_true(length(pima.model.idx) > 0)
  expect_true(length(pima.predict.idx) > 0)
  expect_true(length(pima.model.idx) <= length(pima.predict.idx))
})

testthat::test_that("Test functionality on ames data set, with lots of categorical data", {

  data(ames)
  ames.sample <- 1:(0.75*nrow(ames))
  ames <- prim.data.prepare(ames)
  ames.model <- prim(SalePrice ~ . - PID - Order, data = ames[ames.sample,], peeling.quantile = 0.1, min.support = 0.1)
  ames.predict <- predict(ames.model, ames[-ames.sample,])

  expect_is(ames.model, "prim.peel")
  expect_is(ames.predict, "prim.predict")

  ames.model.idx <- prim.box.index(ames.model, ames)
  ames.predict.idx <- prim.box.index(ames.predict, ames)

  expect_true(length(ames.model.idx) > 0)
  expect_true(length(ames.predict.idx) > 0)
  expect_true(length(ames.model.idx) <= length(ames.predict.idx))
})
