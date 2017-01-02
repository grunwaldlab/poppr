context("Informloci tests")

data("Aeut", package = "poppr")
genos <- c("A/A", "A/B", "A/C", "B/B", "B/C", "C/C")

v  <- sample(genos, 100, replace = TRUE)
w  <- c(rep(genos[2], 99), genos[3])           # found by cutoff
x  <- c(rep(genos[1], 98), genos[3], genos[2]) # found by MAF
y  <- c(rep(genos[1], 99), genos[2])           # found by both
z  <- sample(genos, 100, replace = TRUE)
df <- data.frame(v = v, w = w, x = x, y = y, z = z)

rownames(df) <- .genlab("i", 100)

dat <- df2genind(df, sep = "/")

test_that("informloci finds differentiating samples", {
  expect_message(out <- informloci(dat, MAF = 0), "2 loci found with a cutoff")
  expect_equal(nLoc(out), 3)
  expect_equivalent(locNames(out), c("v", "x", "z"))
})

test_that("informloci finds MAF", {
  expect_message(out <- informloci(dat, cutoff = 0), "2 loci found with MAF")
  expect_equal(nLoc(out), 3)
  expect_equivalent(locNames(out), c("v", "w", "z"))
})

test_that("informloci finds both", {
  expect_message(out <- informloci(dat), "Found 3 uninformative loci")
  expect_equal(nLoc(out), 2)
  expect_equivalent(locNames(out), c("v", "z"))
})

test_that("informloci works for Presence/Absence data", {
  skip_on_cran()
  expect_message(informloci(Aeut), "All sites polymorphic")
  expect_message(informloci(Aeut, MAF = 0.5), "Fewer than 2 loci found informative")
})