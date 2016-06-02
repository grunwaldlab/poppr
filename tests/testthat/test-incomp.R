context("compatibility tests")
data("nancycats", package = "adegenet")

test_that("incomp returns true incomparable samples", {
  skip_on_cran()
  mat <- incomp(nancycats)
  expect_is(mat, "matrix")
  expect_equal(dim(mat), c(0L, 0L))
  mat2 <- incomp(nancycats[pop = c(1, 17), loc = c(1, 4)])
  expect_is(mat2, "matrix")
  expect_true(isSymmetric(mat2))
  expect_equal(dim(mat2), c(15L, 15L))
})