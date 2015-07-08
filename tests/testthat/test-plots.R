context("plotting tests")

data("nancycats", package = "adegenet")
data("Pinf", package = "poppr")
nancy <- popsub(nancycats, c(1, 9))

test_that("info_table plots work", {
	skip_on_cran()
	nancy.miss <- info_table(nancy, plot = TRUE, type = "missing")
	np         <- ggplot2::last_plot()
	nancy.count <- info_table(nancy, plot = TRUE, percent = FALSE, plotlab = FALSE)
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
	expect_output(pp$facet, "facet_null\\(\\)")

	expect_output(np$layers[[1]], "geom_tile")
	expect_output(np$layers[[2]], "geom_hline")
	expect_output(np$layers[[3]], "geom_vline")
	expect_output(np$layers[[4]], "geom_text")

	expect_equal(length(nc$layers), 3)
	expect_output(np$facet, "facet_null\\(\\)")
})

test_that("info_table plots ploidy for diploids", {
	skip_on_cran()
	nancy.ploid <- info_table(nancy, plot = TRUE, type = "ploidy")
	np          <- ggplot2::last_plot()
	expect_equal(unique(np$data$Observed_Ploidy), c(NA, 2))
})

test_that("mlg.table produces barplots", {
	skip_on_cran()
	expect_is(mlg.table(Pinf), "matrix")
	pt <- ggplot2::last_plot()
	expect_is(pt, "ggplot")
	expect_equal(names(pt$data), c("Population", "MLG", "count", "fac"))
	expect_output(pt$facet, "facet_wrap\\(Population\\)")
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
	expect_output(pg$facet, "facet_null\\(\\)")
	expect_output(pg$layers[[1]], "geom_boxplot")
	expect_output(pg$layers[[2]], "geom_hline")
	expect_output(pg$layers[[3]], "geom_text")
})

test_that("ia produces histograms", {
	skip_on_cran()
	res    <- ia(nancy, sample = 20, valuereturn = TRUE, quiet = TRUE)
	iaplot <- ggplot2::last_plot()
	expect_is(res, "ialist")
	expect_output(res, "Index")
	plot(res)
	expect_equivalent(iaplot$data, ggplot2::last_plot()$data)

	pres   <- poppr(nancy, sample = 20, quiet = TRUE, total = FALSE)
	poplot <- ggplot2::last_plot()
	expect_is(pres, "popprtable")
	expect_output(pres, "nancy")


	expect_output(iaplot$layers[[1]], "geom_histogram")
	expect_output(iaplot$layers[[2]], "geom_rug")
	expect_output(iaplot$layers[[3]], "geom_vline")
	expect_output(iaplot$layers[[4]], paste0("geom_text", ".+?", signif(res$index["rbarD"], 3)))
	expect_output(iaplot$layers[[5]], paste0("geom_text", ".+?", signif(res$index["p.rD"], 3)))
	expect_output(iaplot$facet, "facet_null\\(\\)")

	expect_output(poplot$layers[[1]], "geom_histogram")
	expect_output(poplot$layers[[2]], "geom_rug")
	expect_output(poplot$layers[[3]], "geom_vline")
	expect_output(poplot$facet, "facet_wrap\\(population\\)")
})

test_that("pair.ia produces a heatmap", {
	skip_on_cran()
	nan.pair <- pair.ia(nancy, quiet = TRUE)
	pplot    <- ggplot2::last_plot()
	expect_is(nan.pair, c("paria", "matrix"))
	
	plot(nan.pair)
	pplot_lim <- ggplot2::last_plot()

	plot(nan.pair, limits = NULL)
	pplot2 <- ggplot2::last_plot()
	
	expect_equivalent(pplot, pplot2)
	expect_that(pplot2, not(is_equivalent_to(pplot_lim)))

	expect_output(pplot$layers[[1]], "geom_tile")
	expect_output(pplot$layers[[1]], "stat_identity")
	expect_output(pplot$layers[[1]], "position_identity")
	expect_output(pplot$facet, "facet_null\\(\\)")
})

test_that("diversity_ci produces boxplots and correct output", {
	skip_on_cran()
	msg <- "Confidence Intervals have been centered around observed statistic"
	set.seed(999)
	expect_message(datdf <- diversity_ci(Pinf, n = 10, raw = FALSE), msg)
	ciplot <- ggplot2::last_plot()
	expect_is(datdf, "popprtable")
	expect_output(datdf, "\\(\\d\\.\\d{3}, \\d\\.\\d{3})")
	expect_is(ciplot, "ggplot")

	expect_output(ciplot$layers[[1]], "geom_boxplot")
	expect_output(ciplot$layers[[2]], "geom_point")
	expect_output(ciplot$layers[[3]], "geom_errorbar")
	expect_output(ciplot$facet, "facet_wrap\\(Index\\)")

})

test_that("diversity_ci produces boxplots for rarefaction", {
	skip_on_cran()
	msg <- "Samples for rarefaction: 38"
	set.seed(999)
	expect_message(datdf <- diversity_ci(Pinf, n = 10, raw = FALSE, rarefy = TRUE), msg)
	ciplot <- ggplot2::last_plot()
	expect_is(datdf, "popprtable")
	expect_output(datdf, "\\(\\d\\.\\d{3}, \\d\\.\\d{3})")
	expect_output(datdf, "NA")
	expect_is(ciplot, "ggplot")

	expect_equal(length(ciplot$layers), 2)
	expect_output(ciplot$layers[[1]], "geom_boxplot")
	expect_output(ciplot$layers[[2]], "geom_point")
	expect_output(ciplot$facet, "facet_wrap\\(Index\\)")
})


test_that("greycurve produces plots", {
	skip_on_cran()
	expect_output(greycurve(), "")
	expect_output(greycurve(scalebar = TRUE), "")
	expect_output(greycurve(1:100), "")
})