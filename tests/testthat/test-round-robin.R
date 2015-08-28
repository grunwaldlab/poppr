context("Round Robin Tests")

set.seed(69)
x <- matrix(sample(LETTERS[1:3], 15, replace = TRUE), nrow = 5, ncol = 3)
suppressWarnings(x <- df2genind(x, ploidy = 1))
rrx_m <- rrmlg(x)
mlg_truth <- structure(c(3L, 2L, 3L, 1L, 1L, 4L, 3L, 4L, 2L, 1L, 2L, 3L, 2L, 
                         3L, 1L), 
                       .Dim = c(5L, 3L), 
                       .Dimnames = list(NULL, c("L1", "L2", "L3"))
                       )
rrx_f <- rraf(x)
freq_truth <- structure(list(L1 = structure(c(0.333333333333333, 0.666666666666667), 
                                            .Names = c("L1.B", "L1.C")), 
                             L2 = structure(c(0.25, 0.75), 
                                            .Names = c("L2.C", "L2.A")), 
                             L3 = structure(c(0.333333333333333, 0.333333333333333, 
                                              0.333333333333333), 
                                            .Names = c("L3.C", "L3.B", "L3.A"))), 
                        .Names = c("L1", "L2", "L3")
                        )

test_that("rrmlg produces the correct mlgs", {
  expect_equivalent(rrx_m, mlg_truth)
})

test_that("rraf produces correct allele frequencies", {
  skip_on_cran()
  expect_equivalent(rrx_f, freq_truth)
  expect_equivalent(rraf(x, "vector"), unlist(freq_truth))
})

test_that("rraf will correct or not correct allele frequencies when asked", {
  skip_on_cran()
  data(monpop)
  monc  <- rraf(monpop)
  monc_vec <- vapply(monc, sum, numeric(1))
  monnc <- rraf(monpop, correction = FALSE)
  monnc_vec <- vapply(monnc, sum, numeric(1))
  
  expect_more_than(sum(monc_vec), sum(monnc_vec))
  expect_equivalent(monnc_vec, rep(1, nLoc(monpop)))
})

test_that("rraf produces a data frame", {
  skip_on_cran()
  rrx_df <- rraf(x, "data.frame")
  expect_is(rrx_df, "data.frame")
  expect_equivalent(names(rrx_df), c("frequency", "locus", "allele"))
})

test_that("rraf calculates per population", {
  skip_on_cran()
  data(Pram)
  rrx_matrix <- rraf(Pram, by_pop = TRUE, correction = FALSE)
  expect_is(rrx_matrix, "matrix")
  expect_equivalent(dim(rrx_matrix), c(nPop(Pram), ncol(tab(Pram))))
  expect_equivalent(rownames(rrx_matrix), popNames(Pram))
})

test_that("psex and pgen internals produce expected results", {
  skip_on_cran()
  # From Parks and Werth, 1993
  x <- "
 Hk Lap Mdh2 Pgm1 Pgm2 X6Pgd2
54 12 12 12 23 22 11
36 22 22 11 22 33 11
10 23 22 11 33 13 13"
  
  xtab <- read.table(text = x, header = TRUE, row.names = 1)
  xgid <- df2genind(xtab[rep(rownames(xtab), as.integer(rownames(xtab))), ], ncode = 1)
  afreq <- c(Hk.1 = 0.167, Hk.2 = 0.795, Hk.3 = 0.038, 
             Lap.1 = 0.190, Lap.2 = 0.798, Lap.3 = 0.012,
             Mdh2.0 = 0.011, Mdh2.1 = 0.967, Mdh2.2 = 0.022,
             Pgm1.2 = 0.279, Pgm1.3 = 0.529, Pgm1.4 = 0.162, Pgm1.5 = 0.029,
             Pgm2.1 = 0.128, Pgm2.2 = 0.385, Pgm2.3 = 0.487,
             X6Pgd2.1 = 0.526, X6Pgd2.2 = 0.051, X6Pgd2.3 = 0.423)
  freqs   <- afreq[colnames(tab(xgid))]
  pNotGen <- psex(xgid, by_pop = FALSE, freq = freqs, G = 45)
  xpgen   <- exp(rowSums(pgen(xgid, by_pop = FALSE, freq = freqs)))
#   pops <- rep(1L, nInd(xgid))
#   pgenmat <- .Call("get_pgen_matrix_genind", xgid, freqs, pops, 1L, PACKAGE = "poppr")
#   dimnames(pgenmat) <- list(indNames(xgid), locNames(xgid))
#   G <- 45 # Number of MLGs observed
#   xpgen <- exp(rowSums(pgenmat))
#   pNotGen <- (1 - xpgen)^G
#   1 - pNotGen
  res <- matrix(c( unique(xpgen), unique(pNotGen)), ncol = 2)
  
  expected_result <- structure(c(4.14726753733799e-05, 0.00192234266749449, 0.000558567714432373, 
                                 0.00186456862059092, 0.0829457730256321, 0.024829127718338), 
                               .Dim = c(3L,2L))
  expect_equivalent(res, expected_result)
  
})

test_that("psex produces a vector", {
  skip_on_cran()
  data(Pram)
  psexpram <- psex(Pram)
  expect_is(psexpram, "numeric")
  expect_equal(length(psexpram), nInd(Pram))
})

test_that("pgen produces a matrix", {
  skip_on_cran()
  data(Pram)
  pgenlog <- pgen(Pram, by_pop = FALSE)
  expect_is(pgenlog, "matrix")
  expect_equivalent(exp(pgenlog), pgen(Pram, by_pop = FALSE, log = FALSE))
})