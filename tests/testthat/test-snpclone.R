context("snpclone coercion tests")

set.seed(999)
gc <- as.snpclone(glSim(100, 0, n.snp.struc = 1e3, ploidy = 2, parallel = FALSE), parallel = FALSE)
indNames(gc) <- paste("sample", seq(nInd(gc)))
mlg.filter(gc, threads = 1L) <- 0.25


test_that("subsetting a snpclone object retains the MLG definitions", {
  idn <- mlg.id(gc)[[2]]
  idl <- indNames(gc) %in% idn
  idi <- which(idl)
  expect_equal(mll(gc[idn]), mll(gc[idl]))
  expect_equal(mll(gc[idi]), mll(gc[idl]))
})