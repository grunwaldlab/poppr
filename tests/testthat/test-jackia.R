context("Jackknife tests")

data("Pram")
mll(Pram) <- "original"
jhc <- Pram[pop = 3] %>% setPop(~YEAR) %>% popsub("2001")

test_that("jack.ia works out of the box", {
  skip_on_cran()
  x <- jack.ia(jhc)
  expect_is(x, "data.frame")
  expect_equal(nrow(x), 999L)
  expect_equal(ncol(x), 2L)
})

test_that("jack.ia balks at incorrect values of n", {
  skip_on_cran()
  expect_error(jack.ia(jhc, n = nInd(jhc)), "n must be less")
  expect_error(jack.ia(jhc, n = 2), "n must be greater")
})

data("Pinf")

test_that("the mean of the distribution is near the observed value of ia", {
  skip_on_cran()
  rdsa  <- ia(Pinf[pop = 1])[[2]]
  set.seed(999)
  jrdsa <- jack.ia(Pinf[pop = 1])[["rbarD"]]
  testthat::expect_equal(rdsa, mean(jrdsa), tol = 1e-3)
})