library(subgroup.discovery)
context("Patient Rule Induction Method")

test_that("Test functionality on credit data set", {

  data(credit)

  X <- credit[,-6]
  y <- credit$class

  p.train <- prim.default(X = X, y = y, peeling.quantile = 0.1, min.support = 0.4)
  expect_is(p.train, "prim.peel.result")

  rule <- prim.superrule.index(p.train, X)
  expect_is(rule, "logical")
  expect_identical(rule, c(F, F, F, F, F, T, T, T, F, T))

})

test_that("Test functionality on credit data set using the formula interface", {

  data(credit)

  p.train <- prim.formula(class ~ ., data = credit, peeling.quantile = 0.1, min.support = 0.4)
  expect_is(p.train, "prim.peel.result")

  expect_is(p.train$formula, "formula")
  expect_identical(p.train$formula, as.formula(class ~ .))

  rule <- prim.superrule.index(p.train, credit)
  expect_is(rule, "logical")
  expect_identical(rule, c(F, F, F, F, F, T, T, T, F, T))

})

test_that("Test functionality on pima data set", {

  data(pima)

  X <- pima[,-9]
  y <- pima$class

  train <- sample(1:nrow(X), nrow(X) * 0.75)

  p.train <- prim.default(X = X[train,], y = y[train], peeling.quantile = 0.01, min.support = 0.1)
  expect_is(p.train, "prim.peel.result")

  p.test <- prim.test(p.train, X[-train,], y[-train])
  expect_is(p.test, "prim.test.result")


})
