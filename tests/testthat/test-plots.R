context("plotting tests")

data("nancycats", package = "adegenet")
data("Pinf", package = "poppr")
nancy <- popsub(nancycats, c(1, 9))
ggversion <- packageVersion("ggplot2")
oldgg <- package_version("1.0.1")

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

	expect_output(print(pp$layers[[1]]), "geom_tile")

	expect_output(print(np$layers[[1]]), "geom_tile")
	expect_output(print(np$layers[[2]]), "geom_hline")
	expect_output(print(np$layers[[3]]), "geom_vline")
	expect_output(print(np$layers[[4]]), "geom_text")

	expect_equal(length(nc$layers), 3)
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
	expect_output(print(pt$layers), "geom_bar")
})

test_that("genotype_curve produces boxplots", {
	skip_on_cran()
	data("partial_clone", package = "poppr")

	pcg <- genotype_curve(partial_clone, sample = 20, quiet = TRUE)
	expect_is(pcg, "matrix")
	pg <- ggplot2::last_plot()
	expect_is(pg, "ggplot")
	expect_equal(names(pg$data), c("sample", "NumLoci", "MLG"))
	expect_output(print(pg$layers[[1]]), "geom_boxplot")
	expect_output(print(pg$layers[[2]]), "geom_hline")
	expect_output(print(pg$layers[[3]]), "geom_text")
})

test_that("genotype_curve shows output if not quiet", {
	skip_on_cran()
	data("partial_clone", package = "poppr")

	expect_output(genotype_curve(partial_clone, sample = 100), "100%")
})

test_that("genotype_curve can take less than m-1 loci", {
  skip_on_cran()
  nc  <- genotype_curve(nancycats, maxloci = 2, quiet = TRUE)
  nc1 <- genotype_curve(nancycats, maxloci = 1, quiet = TRUE)
  expect_equal(ncol(nc), 2L)
  expect_equal(ncol(nc1), 1L)
  expect_error(genotype_curve(nancycats[loc = 1]), "at least two loci")
})

test_that("ia produces histograms", {
	skip_on_cran()
	res    <- ia(nancy, sample = 20, valuereturn = TRUE, quiet = TRUE)
	iaplot <- ggplot2::last_plot()
	expect_is(res, "ialist")
	expect_output(print(res), "Index")
	plot(res)
	expect_equivalent(iaplot$data, ggplot2::last_plot()$data)

	pres   <- poppr(nancy, sample = 20, quiet = TRUE, total = FALSE)
	poplot <- ggplot2::last_plot()
	expect_is(pres, "popprtable")
	expect_output(print(pres), "nancy")

	if (ggversion <= oldgg){
	  expect_output(print(iaplot$layers[[1]]), "geom_histogram")
	} else {
	  expect_output(print(iaplot$layers[[1]]), "geom_bar")
	}
	
	expect_output(print(iaplot$layers[[2]]), "geom_rug")
	expect_output(print(iaplot$layers[[3]]), "geom_vline")
	expect_output(print(iaplot$layers[[4]]), "geom_text")
	expect_output(print(iaplot$layers[[5]]), "geom_text")

	if (ggversion <= oldgg){
	  expect_output(print(poplot$layers[[1]]), "geom_histogram")
	} else {
	  expect_output(print(poplot$layers[[1]]), "geom_bar")
	}
	
	expect_output(print(poplot$layers[[2]]), "geom_rug")
	expect_output(print(poplot$layers[[3]]), "geom_vline")
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
	
	if (ggversion <= package_version("1.0.1")){
	  expect_equivalent(pplot, pplot2)
	  expect_that(pplot2, not(is_equivalent_to(pplot_lim))) 
	} else {
	  pplotlim     <- pplot$scales$get_scales("fill")$get_limits()
	  pplot2lim    <- pplot2$scales$get_scales("fill")$get_limits()
	  pplot_limlim <- pplot_lim$scales$get_scales("fill")$get_limits()
	  expect_equivalent(pplotlim, pplot2lim)
	  expect_false(identical(pplotlim, pplot_limlim))
	}
	
	

	expect_output(print(pplot$layers[[1]]), "geom_tile")
	expect_output(print(pplot$layers[[1]]), "stat_identity")
	expect_output(print(pplot$layers[[1]]), "position_identity")
})

test_that("diversity_ci produces boxplots and correct output", {
	skip_on_cran()
	msg <- "Confidence Intervals have been centered around observed statistic"
	set.seed(999)
	expect_message(datdf <- diversity_ci(Pinf, n = 10, raw = FALSE), msg)
	ciplot <- ggplot2::last_plot()
	expect_is(datdf, "popprtable")
	expect_output(print(datdf), "\\(\\d\\.\\d{3}, \\d\\.\\d{3})")
	expect_is(ciplot, "ggplot")

	expect_output(print(ciplot$layers[[1]]), "geom_boxplot")
	expect_output(print(ciplot$layers[[2]]), "geom_point")
	expect_output(print(ciplot$layers[[3]]), "geom_errorbar")

})

test_that("diversity_ci produces boxplots for rarefaction", {
	skip_on_cran()
	msg <- "Samples for rarefaction: 38"
	set.seed(999)
	expect_message(datdf <- diversity_ci(Pinf, n = 10, raw = FALSE, rarefy = TRUE), msg)
	ciplot <- ggplot2::last_plot()
	expect_is(datdf, "popprtable")
	expect_output(print(datdf), "\\(\\d\\.\\d{3}, \\d\\.\\d{3})")
	expect_output(print(datdf), "NA")
	expect_is(ciplot, "ggplot")

	expect_equal(length(ciplot$layers), 2)
	expect_output(print(ciplot$layers[[1]]), "geom_boxplot")
	expect_output(print(ciplot$layers[[2]]), "geom_point")
})


test_that("greycurve produces plots", {
	skip_on_cran()
	expect_output(greycurve(), NA)
	expect_output(greycurve(scalebar = TRUE), NA)
	expect_output(greycurve(1:100), NA)
})