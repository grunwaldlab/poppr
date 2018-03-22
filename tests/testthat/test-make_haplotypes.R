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