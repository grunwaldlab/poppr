context("MLG class tests")

## This setup is directly from the man page
set.seed(5000)
(x <- sample(10, 20, replace = TRUE))
(m <- new("MLG", x))

test_that("visible accessor works", {
  skip_on_cran()
  expect_output(show(visible(m)), "original")
  visible(m) <- "contracted"
  expect_output(show(m), "contracted")
})

test_that("MLG2df produces a data frame", {
  skip_on_cran()
  mdf <- MLG2df(m)
  expect_is(mdf, "data.frame")
  expect_identical(mdf, m@mlg)
})

test_that("distname returns a name", {
  skip_on_cran()
  expect_identical(distname(m), "diss.dist")
  distname(m) <- substitute("nei.dist")
  expect_identical(distname(m), "nei.dist")
})

test_that("distargs returns a list", {
  skip_on_cran()
  distname(m) <- substitute("diss.dist")
  expect_identical(distargs(m), list())
  distargs(m) <- list(percent = TRUE)
  expect_identical(distargs(m), list(percent = TRUE))
})

test_that("distalgo returns a specified algorithm", {
  skip_on_cran()
  expect_identical(distalgo(m), "farthest_neighbor")
  distalgo(m) <- "average"
  expect_identical(distalgo(m), "average")
})

test_that("cutoff returns the correct cutoff", {
  skip_on_cran()
  zcut <- setNames(c(0.0, 0.0), c("expanded", "contracted"))
  expect_identical(cutoff(m), zcut)
  cutoff(m)["contracted"] <- 0.2
  tcut <- zcut
  tcut["contracted"] <- 0.2
  expect_identical(cutoff(m), tcut)
})