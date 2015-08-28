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