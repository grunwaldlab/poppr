context("MLG crossplot tests")

data("Pram", package = "poppr")
randpop  <- sample(10, nInd(Pram), replace = TRUE)

test_that("crossplot works without strata", {
  skip_on_cran()
  gnostrat <- mlg.crossplot(Pram, temporal = ~YEAR)
  expect_identical(sort(igraph::V(gnostrat)$name), sort(popNames(Pram)))
})

test_that("crossplot works with strata", {
  skip_on_cran()
  gstrat <- mlg.crossplot(Pram, pop = ~SOURCE, temporal = ~YEAR)
  expect_identical(sort(igraph::V(gstrat)$name), sort(popNames(setPop(Pram, ~SOURCE))))
})

test_that("crossplot works with custom pop", {
  skip_on_cran()
  grand <- mlg.crossplot(Pram, pop = randpop, temporal = ~YEAR)
  expect_equal(sort(as.numeric(igraph::V(grand)$name)), 1:10)
})