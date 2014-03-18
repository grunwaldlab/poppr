context("Missing tests")

data(nancycats, package = "adegenet")

test_that("missingno removes genotypes", {
	expect_that(genmiss <- missingno(nancycats, "geno", cutoff = 0.05), prints_text("38 genotypes contained missing values greater than 5%."))
	expect_that(nInd(genmiss), equals(199))
})

test_that("missingno removes loci", {
	expect_that(locmiss <- missingno(nancycats, "loci", cutoff = 0.05), prints_text("2 loci contained missing values greater than 5%."))
	expect_that(nLoc(locmiss), equals(7))
})

test_that("missingno matches na.replace", {
	expect_identical(missingno(nancycats, "zero", quiet = TRUE), na.replace(nancycats, "0", quiet = TRUE))
	expect_identical(missingno(nancycats, "mean", quiet = TRUE), na.replace(nancycats, "mean", quiet = TRUE))
})