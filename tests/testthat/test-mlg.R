context("Multilocus genotype tests")

test_that("multilocus genotype vector is same length as samples", {
  data(Aeut, package = "poppr")
  data(partial_clone, package = "poppr")
  data(nancycats, package = "adegenet")
  amlg <- mlg.vector(Aeut)
  pmlg <- mlg.vector(partial_clone)
  nmlg <- mlg.vector(nancycats)
  expect_that(length(amlg), equals(nInd(Aeut)))
  expect_that(length(pmlg), equals(nInd(partial_clone)))
  expect_that(length(nmlg), equals(nInd(nancycats)))
  expect_that(length(unique(amlg)), equals(mlg(Aeut, quiet = TRUE)))
  expect_that(length(unique(pmlg)), equals(mlg(partial_clone, quiet = TRUE)))
  expect_that(length(unique(nmlg)), equals(mlg(nancycats, quiet = TRUE)))
})

test_that("multilocus genotype matrix matches mlg.vector and data", {
  data(Aeut, package = "poppr")
  data(partial_clone, package = "poppr")
  data(nancycats, package = "adegenet")
  aclone <- as.genclone(Aeut)
  atab   <- mlg.table(Aeut, bar = FALSE)
  ptab   <- mlg.table(partial_clone, bar = FALSE)
  ntab   <- mlg.table(nancycats, bar = FALSE)
  expect_that(nrow(atab), equals(length(Aeut@pop.names)))
  expect_that(nrow(ptab), equals(length(partial_clone@pop.names)))
  expect_that(nrow(ntab), equals(length(nancycats@pop.names)))
  expect_that(ncol(atab), equals(mlg(Aeut, quiet = TRUE)))
  expect_that(ncol(ptab), equals(mlg(partial_clone, quiet = TRUE)))
  expect_that(ncol(ntab), equals(mlg(nancycats, quiet = TRUE)))
  expect_that(sum(atab), equals(nInd(Aeut)))
  expect_that(sum(ptab), equals(nInd(partial_clone)))
  expect_that(sum(ntab), equals(nInd(nancycats)))
})

test_that("mlg.crosspop will work with subsetted genclone objects", {
  data(Aeut, package = "poppr")
  agc <- as.genclone(Aeut)
  Athena <- popsub(agc, "Athena")
  setpop(Athena) <- ~Subpop
  expected_output <- structure(list(MLG.13 = structure(c(1L, 1L), .Names = c("8", 
"9")), MLG.23 = structure(c(1L, 1L), .Names = c("4", "6")), MLG.24 = structure(c(1L, 
1L), .Names = c("9", "10")), MLG.32 = structure(c(1L, 1L), .Names = c("7", 
"9")), MLG.52 = structure(c(1L, 1L), .Names = c("5", "9")), MLG.63 = structure(c(1L, 
1L), .Names = c("1", "5"))), .Names = c("MLG.13", "MLG.23", "MLG.24", 
"MLG.32", "MLG.52", "MLG.63"))
  expected_mlgout <- c(13, 23, 24, 32, 52, 63)

  expect_that(x <- mlg.crosspop(Athena, quiet = TRUE), equals(expected_output))
  expect_that(y <- mlg.crosspop(Athena, indexreturn = TRUE), equals(expected_mlgout))
})




