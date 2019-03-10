context("Resampling tests")

data("Pram")
mll(Pram) <- "original"
jhc <- Pram[pop = 3] %>% setPop(~YEAR) %>% popsub("2001")

test_that("resample.ia works out of the box", {
  skip_on_cran()
  x <- resample.ia(jhc)
  expect_is(x, "data.frame")
  expect_equal(nrow(x), 999L)
  expect_equal(ncol(x), 2L)
})

test_that("boot.ia works out of the box", {
  skip_on_cran()
  x <- boot.ia(jhc)
  expect_is(x, "data.frame")
  expect_equal(nrow(x), 999L)
  expect_equal(ncol(x), 2L)
})

test_that("resample and boot balk at incorrect values of n", {
  skip_on_cran()
  expect_error(resample.ia(jhc, n = nInd(jhc)), "n must be fewer")
  expect_error(resample.ia(jhc, n = 2), "n must be greater")
  expect_error(boot.ia(jhc[1:2], how = "full"), "n must be greater")
})

data("Pinf")

test_that("the mean of the distribution is near the observed value of ia", {
  skip_on_cran()
  rdsa  <- ia(Pinf[pop = 1])[[2]]
  set.seed(999)
  jrdsa <- resample.ia(Pinf[pop = 1])[["rbarD"]]
  testthat::expect_equal(rdsa, mean(jrdsa), tol = 1e-2)
})

test_that("jack.ia is deprecated", {
  skip_on_cran()
  expect_warning(x <- jack.ia(Pinf, reps = 9, quiet = TRUE), "jack.ia\\(\\) is deprecated")
  expect_is(x, "data.frame")
})
