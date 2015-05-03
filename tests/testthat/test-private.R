context("Private Allele Tests")

data("Pinf", package = "poppr")

test_that("Private alleles can be calculated for hierarchies", {
  by_continent <- private_alleles(Pinf, all ~ Continent)
  by_country   <- private_alleles(Pinf, all ~ Country)
  expect_true(nrow(by_continent) <= nPop(setPop(Pinf, ~Continent)))
  expect_true(nrow(by_country) <= nPop(setPop(Pinf, ~Country)))
})

test_that("Private alleles can be tabluated over loci", {
  locus_dose <- private_alleles(Pinf, locus ~ .)
  locus_nodose <- private_alleles(Pinf, locus ~ ., count.alleles = FALSE)
  expect_that(ncol(locus_dose), equals(ncol(locus_nodose)))
  expect_true(all(locus_nodose <= locus_dose))
})