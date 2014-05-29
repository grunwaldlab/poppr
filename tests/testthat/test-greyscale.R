context("Greyscale tests")

test_that("rerange works as expected", {
	xnorm   <- rnorm(1000)
	xrange  <- seq(0, 1, length = 100)
	xrange2 <- seq(0, 2, length = 100)
	expect_that(range(poppr:::rerange(xnorm)), equals(c(0, 1)))
	expect_that(range(poppr:::rerange(xrange)), equals(c(0, 1)))
	expect_that(range(poppr:::rerange(xrange2)), equals(c(0, 1)))
})