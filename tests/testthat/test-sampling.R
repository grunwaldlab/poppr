context("Random Sample Tests")

data("nancycats", package = "adegenet")
data("Pinf", package = "poppr")

test_that("all shuffling methods work with missing data", {
	skip_on_cran()
	nan1   <- popsub(nancycats, 1)
	inftab <- info_table(nan1)

	expect_error(s1 <- shufflepop(nan1, method = 1), NA)
	expect_error(s2 <- shufflepop(nan1, method = 2), NA)
	expect_error(s3 <- shufflepop(nan1, method = 3), NA)
	expect_error(s4 <- shufflepop(nan1, method = 4), NA)

	expect_that(info_table(s1), is_equivalent_to(inftab))
	expect_null(suppressMessages(info_table(s2)))
	expect_null(suppressMessages(info_table(s3)))
	expect_that(info_table(s4), is_equivalent_to(inftab))
})

test_that("shuffling methods behave as expected with polyploids", {
	skip_on_cran()
	pr <- recode_polyploids(Pinf, newploidy = TRUE)

	expect_error(s1 <- shufflepop(pr, method = 1), NA)
	expect_error(s2 <- shufflepop(pr, method = 2), NA)
	expect_error(s3 <- shufflepop(pr, method = 3), NA)
	expect_error(s4 <- shufflepop(pr, method = 4), NA)

	s2rows <- rowSums(tab(s2), na.rm = TRUE)
	s3rows <- rowSums(tab(s3), na.rm = TRUE)
	expect_that(rowSums(tab(s1)), equals(rowSums(tab(pr))))
	expect_true(all(s2rows == sample(s2rows, 1)))
	expect_true(all(s3rows == sample(s3rows, 1)))
	expect_that(rowSums(tab(s4)), not(equals(rowSums(tab(pr)))))

})

test_that("shuffling methods work for PA data", {
	skip_on_cran()
	data("Aeut", package = "poppr")
	expect_is(shufflepop(Aeut, method = 1), "genind")
	expect_is(shufflepop(Aeut, method = 2), "genind")
	expect_is(shufflepop(Aeut, method = 3), "genind")
	expect_is(shufflepop(Aeut, method = 4), "genind")

	A10 <- Aeut[sample(nInd(Aeut), 10)]
	expect_is(poppr(A10, sample = 9, method = 1, quiet = TRUE, hist = FALSE, sublist = "Total"), "popprtable")
	expect_is(poppr(A10, sample = 9, method = 2, quiet = TRUE, hist = FALSE, sublist = "Total"), "popprtable")
	expect_is(poppr(A10, sample = 9, method = 3, quiet = TRUE, hist = FALSE, sublist = "Total"), "popprtable")
	expect_is(poppr(A10, sample = 9, method = 4, quiet = TRUE, sublist = "Total"), "popprtable")

})
