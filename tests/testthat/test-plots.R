context("plotting tests")

data("nancycats", package = "adegenet")
data("Pinf", package = "poppr")

test_that("info_table plots work", {
	skip_on_cran()
	nancy.miss <- info_table(nancycats, plot = TRUE, type = "missing")
	np         <- ggplot2::last_plot()
	nancy.count <- info_table(nancycats, plot = TRUE, percent = FALSE, plotlab = FALSE)
	nc          <- ggplot2::last_plot()
	Pinf.ploid <- info_table(Pinf, plot = TRUE, type = "ploidy")
	pp         <- ggplot2::last_plot()
	expect_is(nancy.miss, "locustable")
	expect_equal(names(np$data), c("Locus", "Population", "Missing"))
	expect_is(Pinf.ploid, "locustable")
	expect_equal(names(pp$data), c("Samples", "Loci", "Observed_Ploidy"))
	expect_equivalent(nancy.miss[1, 1], 2/10)
	expect_equivalent(nancy.count[1, 1], 2)
	expect_equal(Pinf.ploid["PiCO01", "Pi4B"], 3)

	expect_output(pp$layers[[1]], "geom_tile")
	expect_output(pp$facet, "facet_null()")

	expect_output(np$layers[[1]], "geom_tile")
	expect_output(np$layers[[2]], "geom_hline")
	expect_output(np$layers[[3]], "geom_vline")
	expect_output(np$layers[[4]], "geom_text")

	expect_equal(length(nc$layers), 3)
	expect_output(np$facet, "facet_null()")
})

test_that("info_table plots ploidy for diploids", {
	skip_on_cran()
	nancy.ploid <- info_table(nancycats, plot = TRUE, type = "ploidy")
	np          <- ggplot2::last_plot()
	expect_equal(unique(np$data$Observed_Ploidy), c(NA, 2))
})

test_that("mlg.table produces barplots", {
	skip_on_cran()
	expect_is(mlg.table(Pinf), "matrix")
	pt <- ggplot2::last_plot()
	expect_is(pt, "ggplot")
	expect_equal(names(pt$data), c("Population", "MLG", "count", "fac"))
	expect_output(pt$facet, "facet_wrap.Population.")
	expect_output(pt$layers, "geom_bar")
})

test_that("genotype_curve produces boxplots", {
	skip_on_cran()
	data("partial_clone", package = "poppr")

	pcg <- genotype_curve(partial_clone, sample = 20, quiet = TRUE)
	expect_is(pcg, "matrix")
	pg <- ggplot2::last_plot()
	expect_is(pg, "ggplot")
	expect_equal(names(pg$data), c("sample", "NumLoci", "MLG"))
	expect_output(pg$facet, "facet_null()")
	expect_output(pg$layers[[1]], "geom_boxplot")
	expect_output(pg$layers[[2]], "geom_hline")
	expect_output(pg$layers[[3]], "geom_text")
})

test_that("ia produces histograms", {
	skip_on_cran()
	nan5 <- popsub(nancycats, 5)
	res <- ia(nan5, sample = 99, valuereturn = TRUE, quiet = TRUE)
	iaplot <- ggplot2::last_plot()
	expect_is(res, "ialist")
	plot(res)
	expect_equivalent(iaplot$data, ggplot2::last_plot()$data)
	
	expect_output(iaplot$layers[[1]], "geom_histogram")
	expect_output(iaplot$layers[[2]], "geom_rug")
	expect_output(iaplot$layers[[3]], "geom_vline")
	expect_output(iaplot$layers[[4]], paste0("geom_text", ".+?", signif(res$index["rbarD"], 3)))
	expect_output(iaplot$layers[[5]], paste0("geom_text", ".+?", signif(res$index["p.rD"], 3)))
})