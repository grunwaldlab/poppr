context("Polyploid Tests")

test_that("recode_polyploids works as expected", {
	testdf  <- data.frame(test = c("00/20/23/24", "20/24/26/43"))
	testgid <- df2genind(testdf, ploidy = 4, sep = "/")
	testrec <- recode_polyploids(testgid)
	expect_equal(testgid@ploidy, testrec@ploidy)
	expect_equivalent(testgid@tab[2, -1], testrec@tab[2, ])
	expect_equivalent(testrec@tab[1, ], c(1/3, 1/3, 1/3, 0, 0))
})

test_that("recode_polyploids won't take diploid data", {
	data(nancycats, package = "adegenet")
	expect_warning(recode_polyploids(nancycats), "Input is not a polyploid data set, returning original.")
})