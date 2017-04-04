context("snpclone coercion tests")

set.seed(999)
dat            <- matrix(sample(c(0:2), 50*64, replace = TRUE), 50, 64)
dat[sample(length(dat), 10)] <- NA
dat2           <- rbind(dat, tail(dat, 10))
rownames(dat2) <- paste("sample", seq(60))
colnames(dat2) <- paste("SNP", seq(64))
gl             <- new("genlight", dat2, parallel = FALSE, n.cores = 1L)
sc             <- as.snpclone(gl, parallel = FALSE, n.cores = 1L)
# mlg.filter(sc, threads = 1L) <- 0.25


test_that("subsetting a snpclone object retains the MLG definitions", {
  idn <- mlg.id(sc)[[2]]
  idl <- indNames(sc) %in% idn
  idi <- which(idl)
  expect_equal(mll(sc[idn]), mll(sc[idl]))
  expect_equal(mll(sc[idi]), mll(sc[idl]))
})

test_that("mlg.vector returns the same value for genlight and snpclone", {
  skip_on_cran()
  expect_equal(mlg.vector(gl), mlg.vector(sc))
})

test_that("snpclone objects can be initialized with unique mlgs", {
  skip_on_cran()
  scu <- as.snpclone(gl, parallel = FALSE, n.cores = 1L, mlg = seq_len(nInd(gl)))
  expect_identical(indNames(scu), indNames(sc))
  expect_identical(locNames(scu), locNames(sc))
  expect_identical(mll(scu), seq_len(nInd(gl)))
  expect_identical(mll(scu[mlg.reset = TRUE]), mll(sc))
})
# TODO: test for handling of missing data (should be considered new genotype)
test_that("a single instance of missing data won't induce a new genotype", {
  skip_on_cran()
  datm        <- dat2
  datm[60, 1] <- NA
  scm         <- as.snpclone(new("genlight", datm, parallel = FALSE, n.cores = 1L),
                             parallel = FALSE, n.cores = 1L)
  expect_equal(nmll(scm), 50L)
})

test_that("mlg.filter will reduce the number of multilocus lineages", {
  skip_on_cran()
  mlg.filter(sc) <- 0.4
  expect_equal(nmll(sc), 24)
})

test_that("mlg.filter<- can't be used for genlight objects", {
  skip_on_cran()
  expect_warning(mlg.filter(gl) <- 0.5, "mlg.filter<- only has an effect on snpclone objects.")
})