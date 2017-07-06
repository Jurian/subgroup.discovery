library(subgroup.discovery)
context("Patient Rule Induction Method")

testthat::test_that("Test functionality on credit data set", {

  data(credit)

  X <- credit[,-6]
  y <- credit$class

  p.train <- subgroup.discovery::prim.default(X = X, y = y, peeling.quantile = 0.1, min.support = 0.4)
  expect_is(p.train, "prim.peel")

  rule <- subgroup.discovery::prim.superrule.index(p.train, X)
  expect_is(rule, "logical")
  expect_identical(rule, c(F, F, F, F, F, T, T, T, F, T))

})

testthat::test_that("Test functionality on credit data set using the formula interface", {

  data(credit)

  p.train <- subgroup.discovery::prim.formula(class ~ ., data = credit, peeling.quantile = 0.1, min.support = 0.4)
  expect_is(p.train, "prim.peel")

  expect_is(p.train$formula, "formula")
  expect_identical(p.train$formula, as.formula(class ~ .))

  rule <- subgroup.discovery::prim.superrule.index(p.train, credit)
  expect_is(rule, "logical")
  expect_identical(rule, c(F, F, F, F, F, T, T, T, F, T))

})

testthat::test_that("Test the validating functionality on pima data set", {

  data(pima)

  X <- pima[,-9]
  y <- pima$class

  train <- sample(1:nrow(X), nrow(X) * 0.66)

  p.train <- subgroup.discovery::prim.default(X = X[train,], y = y[train], peeling.quantile = 0.05, min.support = 0.1)
  expect_is(p.train, "prim.peel")

  p.test <- subgroup.discovery::prim.validate(p.train, X[-train,], y[-train])
  expect_is(p.test, "prim.validate")

})

testthat::test_that("Test the PRIM covering algorithm using the pima data set", {

  data(pima)

  X <- pima[,-9]
  y <- pima$class

  p.cov <- subgroup.discovery::prim.cover.default(X, y, peeling.quantile = 0.05, min.support = 0.1, max.boxes = 3)

  expect_is(p.cov, "prim.cover")
  #expect_identical(length(p.cov$covers), 3)
  expect_true(!is.unsorted(rev(sapply(p.cov$covers, function(x) x$cov.N))))

})

testthat::test_that("Test the PRIM covering algorithm using the formula interface", {

  data(pima)

  p.cov <- subgroup.discovery::prim.cover.formula(class ~ ., pima, peeling.quantile = 0.05, min.support = 0.1, max.boxes = 3)

  expect_is(p.cov, "prim.cover")
  #expect_identical(length(p.cov$covers), 3)
  expect_true(!is.unsorted(rev(sapply(p.cov$covers, function(x) x$cov.N))))

})
