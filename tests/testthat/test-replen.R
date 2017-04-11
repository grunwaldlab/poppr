context("repeat length handling")

data("nancycats", package = "adegenet")
names(alleles(nancycats)) <- locNames(nancycats)
nanrep       <- rep(2, 9)
named_nanrep <- setNames(nanrep, locNames(nancycats))
nantest       <- test_replen(nancycats, nanrep)
nanfix        <- fix_replen(nancycats, nanrep)
named_nantest <- test_replen(nancycats, named_nanrep)
named_nanfix  <- fix_replen(nancycats, named_nanrep)
test_that("test_replen and fix_replen works as expected for conguent vectors", {
  expect_equal(sum(nantest), 5)
  expect_identical(nantest, named_nantest)
  expect_identical(nanfix, named_nanfix)
  expect_true(all(floor(nanfix)[!nantest] == 1))
})

test_that("test_replen and fix_replen work for larger length vectors", {
  nanrep10       <- c(nanrep, 5)
  named_nanrep11 <- c(named_nanrep, foo = 5, bar = 5)
  bad_named_nanrep11 <- named_nanrep11
  names(bad_named_nanrep11)[1] <- "bob"
  expect_error(test_replen(nancycats, nanrep10), "length of repeats \\(10\\)")
  expect_error(fix_replen(nancycats, nanrep10), "length of repeats \\(10\\)")
  expect_warning(expect_error(test_replen(nancycats, bad_named_nanrep11), "repeat lengths... bob, foo, bar"), "bob")
  expect_warning(nts <- test_replen(nancycats, named_nanrep11), "foo, bar")
  expect_warning(nfx <- fix_replen(nancycats, named_nanrep11), "foo, bar")
  expect_identical(nts, nantest)
  expect_identical(nfx, nanfix)
})

test_that("test_replen and fix_replen will not work for short vectors", {
  expect_error(test_replen(nancycats, nanrep[1:7]), "length of repeats \\(7\\)")
  expect_error(test_replen(nancycats, named_nanrep[1:7]), "length of repeats \\(7\\)")
})

test_that("fix_replen will work for shuffled vectors", {
  data(Pram)
  (Pram_replen <- setNames(c(3, 2, 4, 4, 4), locNames(Pram)))
  shuff <- fix_replen(Pram, sample(Pram_replen))
  orig  <- fix_replen(Pram, Pram_replen)
  expect_equal(shuff, orig)
})

test_that("fix_replen throws errors for weird replens", {
  skip_on_cran()
  data(partial_clone)
  expect_warning(fix_replen(partial_clone, rep(10, 10)), paste(locNames(partial_clone), collapse = ", "))
  expect_warning(fix_replen(partial_clone, rep(10, 10), fix_some = FALSE), "Original repeat lengths are being returned")
  expect_warning(fix_replen(partial_clone, rep(2, 10)), "The repeat lengths for Locus_2, Locus_7, Locus_9 are not consistent.")
  expect_warning(fix_replen(partial_clone, rep(2, 10)), "Repeat lengths with some modification are being returned: Locus_3")
})
