context("snpclone coercion tests")

set.seed(999)
gl <- glSim(50, 0, n.snp.struc = 1e3, ploidy = 2, parallel = FALSE)

indNames(gl) <- paste("sample", seq(nInd(gl)))
sc           <- as.snpclone(gl, parallel = FALSE)
mlg.filter(sc, threads = 1L) <- 0.25


test_that("subsetting a snpclone object retains the MLG definitions", {
  idn <- mlg.id(sc)[[2]]
  idl <- indNames(sc) %in% idn
  idi <- which(idl)
  expect_equal(mll(sc[idn]), mll(sc[idl]))
  expect_equal(mll(sc[idi]), mll(sc[idl]))
})

# TODO: test mlg.vector on both gl and sc. They should be the same

# TODO: test for mlg.reset argument by assigning mlg = seq_len(nInd(gl)) in as.snpclone

# TODO: test for handling of missing data (should be considered new genotype)

# TODO: test to see if missing data at the same locus is considered similar or different