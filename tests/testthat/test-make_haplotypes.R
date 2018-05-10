context("Haplotype Splitting")

data(partial_clone)

test_that("make_haplotypes returns double the number of samples for diploids", {
  expect_warning(pc <- make_haplotypes(partial_clone[seq(5)]), "strata")
  expect_equal(nInd(pc), 10)
})

test_that("make_haplotypes returns haploids unharmed", {
  data(monpop)
  expect_warning(mon <- make_haplotypes(monpop), "haploid data")
  expect_identical(mon, monpop)
})

test_that("make_haplotypes only takes genind objects", {
  expect_error(make_haplotypes(1:10), "genind")
})

test_that("make_haplotypes will work for genlight objects", {
  skip_on_cran()
  set.seed(999)
  gl <- glSim(10, 10, 10, ploidy = 2)
  gl1 <- glSim(10, 10, 10, ploidy = 1)
  # Strata are added automatically
  expect_warning(res <- make_haplotypes(gl), "No strata")
  expect_equal(nInd(res), 20L)
  strata(gl) <- data.frame(pop = pop(gl))
  res2       <- make_haplotypes(gl)
  expect_identical(res, res2)
  # haploids are rejected
  expect_warning(hap <- make_haplotypes(gl1), "haploid")
  expect_identical(hap, gl1)
})