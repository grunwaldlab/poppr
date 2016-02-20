context("Errors and Exceptions")

data("partial_clone", package = "poppr")
reps <- rep(1, 10)

test_that("bruvo.msn throws errors with incorrect add and loss arguments", {
  expect_error(bruvo.msn(partial_clone[1:10], replen = reps, add = "1", loss = TRUE))
  expect_error(bruvo.msn(partial_clone[1:10], replen = reps, add = TRUE, loss = "TRUE"))
  expect_error(bruvo.msn(partial_clone[1:10], replen = reps, add = matrix(TRUE, 10, 10), loss = TRUE))
})

test_that("bruvo.dist throws an error when given non-microsatellite data", {
	data("Aeut", package = "poppr")
	data("H3N2", package = "adegenet")
	expect_error(bruvo.dist(Aeut))
	expect_error(bruvo.dist(H3N2))
})

test_that("bruvo.dist estimates repeat length when replen is not present", {
	expect_warning(bruvo.dist(partial_clone), "c\\(3, 1, 2, 1, 2, 4, 2, 3, 1, 2\\)")
})

test_that("poppr rejects genlight objects", {
  skip_on_cran()
  dat <- list(toto=c(1,1,0,0), titi=c(NA,1,1,0), tata=c(NA,0,3, NA))
  x   <- new("genlight", dat, parallel = FALSE)
  expect_error(poppr(x))
})