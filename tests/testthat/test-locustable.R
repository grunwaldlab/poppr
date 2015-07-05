context("Locus table tests")

data("nancycats", package = "adegenet")
data("Pinf", package = "poppr")
nancy <- popsub(nancycats, c(1, 9))

test_that("locus_table correctly treats polyploids", {
	skip_on_cran()
	pinflt <- locus_table(Pinf)
	pinfpp <- poppr(Pinf, quiet = TRUE)
	expect_equivalent(pinflt[-(nLoc(Pinf) + 1), "allele"], nAll(Pinf) - 1)
	expect_equivalent(pinflt["mean", "Hexp"], pinfpp[pinfpp$Pop == "Total", "Hexp"])
	salt   <- locus_table(Pinf, population = "South America")
	expect_equivalent(salt["Pi33", "allele"], 1)
	expect_equivalent(salt["mean", "Hexp"], pinfpp[pinfpp$Pop == "South America", "Hexp"])
})

test_that("locus_table presents different stats", {
	skip_on_cran()
	expect_message(nanlt <- locus_table(nancy), "Simpson")
	expect_message(nanlt <- locus_table(nancy, index = "shannon"), "Shannon")
	expect_message(nanlt <- locus_table(nancy, index = "invsimpson"), "Taylor")
	expect_message(nangt <- locus_table(nancy, lev = "genotype"), "genotype")
	expect_more_than(nangt[, "genotype"], nanlt[, "allele"])
	expect_output(nanlt <- locus_table(nancy, information = FALSE), "")
})