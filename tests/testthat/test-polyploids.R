context("Polyploid Tests")

test_that("recode_polyploids works as expected", {
	testdf  <- data.frame(test = c("00/20/23/24", "20/24/26/43"))
	testgid <- df2genind(testdf, ploidy = 4, sep = "/")
	testrec <- recode_polyploids(testgid)
	expect_equal(testgid@ploidy, testrec@ploidy)
	expect_equivalent(testgid@tab[2, -1], testrec@tab[2, ])
	expect_equivalent(testrec@tab[1, ], c(1, 1, 1, 0, 0))
})

test_that("recode_polyploids won't take diploid data", {
	data(nancycats, package = "adegenet")
	expect_warning(recode_polyploids(nancycats), "Input is not a polyploid data set, returning original.")
})

test_that("recode_polyploids will go there and back again", {
	skip_on_cran()
	data(Pinf, package = "poppr")
	expect_warning(pr <- recode_polyploids(Pinf, addzero = TRUE), "addzero = TRUE")
	pr <- recode_polyploids(Pinf, newploidy = TRUE)
	expect_output(show(pr), "86 diploid \\(55\\) and triploid \\(31\\) individuals")
	expect_that(length(unique(ploidy(pr))), equals(2))
	zpr <- recode_polyploids(pr, addzero = TRUE)
	expect_that(length(unique(ploidy(zpr))), equals(1))
	przpr <- recode_polyploids(zpr, newploidy = TRUE)
	rectab <- tab(przpr)
	expect_equivalent(rectab[, colnames(tab(pr))], tab(pr))
})


test_that("as.genambig will convert correctly", {
  skip_on_cran()
  testthat::skip_if_not_installed("polysat")
  data(Pinf, package = "poppr")
  expect_is(ps <- as.genambig(Pinf), "genambig")
})

test_that("haplodiploids will work", {
  skip_on_cran()
  df <- matrix(c("0/2", "2/1", "3/0", 
                 "0/1", "2/2", "1/1"), nrow = 3)
  x <- df2genind(df, sep = "/")
  expect_equal(ploidy(x), rep(2L, 3))
  y <- recode_polyploids(x, newploidy = TRUE)
  expect_equal(ploidy(y), c(1L, 2L, 2L))
  z <- recode_polyploids(y, addzero = TRUE)
  expect_equal(tab(x), tab(z))
})