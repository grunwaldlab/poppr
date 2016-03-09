context("Missing tests")


mat <- matrix(sample(5, 100, replace = TRUE), 10, 10)
mat[lower.tri(mat)] <- NA
dat <- adegenet::df2genind(mat, ploidy = 1)

test_that("missingno removes genotypes", {
  expect_message(genmiss <- missingno(dat, "geno", cutoff = 0.5), "4 genotypes contained missing values greater than 50%")
  expect_that(nInd(genmiss), equals(6))
})

test_that("missingno removes loci", {
  expect_message(locmiss <- missingno(dat, "loci", cutoff = 0.5), "4 loci contained missing values greater than 50%")
  expect_that(nLoc(locmiss), equals(6))
})

test_that("missingno matches na.replace", {
  mnz <- tab(missingno(dat, "zero", quiet = TRUE))
  mnm <- tab(missingno(dat, "mean", quiet = TRUE))
  expect_identical(mnz, tab(dat, NA.method = "zero", quiet = TRUE))
  expect_identical(mnm, tab(dat, NA.method = "mean", quiet = TRUE))
})