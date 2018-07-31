# Setup -------------------------------------------------------------------

data("nancycats", package = "adegenet")
data("Pinf", package = "poppr")
data("Pram", package = "poppr")
nancy <- popsub(nancycats, c(1, 9))
ggversion <- packageVersion("ggplot2")
oldgg <- package_version("1.0.1")
pcap <- function(x) print(capture.output(x))
mll(Pram) <- "original"

# info_table plots --------------------------------------------------------


context("info_table plots")

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

test_that("info_table also works on non-SSR data", {
  skip_on_cran()
  data("H3N2", package = "adegenet")
  expect_failure(expect_error(info_table(H3N2[50:60, loc = 1:10], type = "ploidy")))
})

context("mlg.table plots")

# mlg.table plots ---------------------------------------------------------


test_that("mlg.table produces barplots", {
	skip_on_cran()
  pm <- mlg.table(Pinf)
  pt <- ggplot2::last_plot()
	expect_is(pm, "matrix")
	expect_identical(pm, mlg.table(Pinf, color = TRUE))
	ptc <- ggplot2::last_plot()
	expect_identical(pm, mlg.table(Pinf, background = TRUE))
	ptb <- ggplot2::last_plot()
	expect_is(pt, "ggplot")
	expect_is(ptc, "ggplot")
	expect_is(ptb, "ggplot")
	expect_equal(names(pt$data), c("Population", "MLG", "count", "order"))
	expect_equal(names(ptc$data), c("Population", "MLG", "count", "order", "n"))
	expect_equal(names(ptb$data), c("Population", "MLG", "count", "order", "n"))
	
	expect_identical(ptb$data, ptc$data)
	expect_false(identical(ptb$data, pt$data))
	
	expect_output(print(pt$layers), "geom_bar")
})

test_that("mlg.table will plot color plot without total", {
  skip_on_cran()
  invisible(mlg.table(Pinf, total = TRUE, color = TRUE))
  pttc <- ggplot2::last_plot()
  p    <- unique(pttc$data$Population) %>% 
    as.character() %>% 
    sort()
  expect_equal(p, sort(popNames(Pinf)))
})

test_that("mlg.table will utilize old versions of dplyr", {
  skip_on_cran()
  options(poppr.old.dplyr = TRUE)
  expect_silent(x <- mlg.table(Pinf, background = TRUE))
  expect_silent(x <- mlg.table(Pinf))
  options(poppr.old.dplyr = FALSE)
})


context("genotype_curve plots")

# genotype_curve plots ----------------------------------------------------



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
	pd <- getOption("poppr.debug")
	options(poppr.debug = TRUE)
	expect_output(genotype_curve(partial_clone, sample = 100), "100%")
	options(poppr.debug = pd)
})

test_that("genotype_curve can take less than m-1 loci", {
  skip_on_cran()
  nc  <- genotype_curve(nancycats, maxloci = 2, quiet = TRUE)
  nc1 <- genotype_curve(nancycats, maxloci = 1, quiet = TRUE)
  expect_equal(ncol(nc), 2L)
  expect_equal(ncol(nc1), 1L)
  expect_error(genotype_curve(nancycats[loc = 1]), "at least two loci")
})

test_that("genotype_curve will not plot if you tell it", {
  skip_on_cran()
  # No output here
  expect_output(pcap(genotype_curve(nancycats, maxloci = 1, quiet = TRUE)), "character\\(0\\)")
  # Some output here
  expect_output(pcap(genotype_curve(nancycats, maxloci = 1, quiet = TRUE, plot = FALSE)))
})

context("ia plots")

# ia plots ----------------------------------------------------------------


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

	if (ggversion <= oldgg) {
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


test_that("ia will still plot if the observed value is NA or NaN", {
  skip_on_cran()
  res <- ia(Pram[pop = 9], sample = 99, quiet = TRUE)
  expect_true(is.nan(res["rbarD"]))
  expect_true(is.na(res["p.rD"]))
})


test_that("ia will know where to place the label depending on the mean", {
  skip_on_cran()
  set.seed(99)
  res <- ia(Pram[pop = 7], sample = 99, quiet = TRUE)
  p <- ggplot2::last_plot()
  prange <- range(p$data$value)
  ptext  <- p$layers[[4]]$data$x
  expect_true(ptext %in% prange)
  expect_false(identical(res[["rbarD"]], ptext))
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

test_that("pair.ia can plot p-values", {
  skip_on_cran()
  p1 <- ggplot2::last_plot()
  tmp <- matrix(round(runif(45*4), 3), nrow = 45, ncol = 4)
  rownames(tmp) <- apply(combn(letters[1:10], 2), 2, paste, collapse = ":")
  colnames(tmp) <- c("Ia", "p.Ia", "rbarD", "p.rD")
  class(tmp) <- c("pairia", "matrix")
  plot(tmp)
  p2 <- ggplot2::last_plot()
  expect_is(p1, "ggplot")
  expect_is(p2, "ggplot")
  # The previous plot is not the same as the current plot
  expect_failure(expect_identical(p1, p2))
})

context("diversity_ci plots")

# diversity_ci plots ------------------------------------------------------


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

context("greycurve plots")

# greycurve plots ---------------------------------------------------------


test_that("greycurve produces plots", {
	skip_on_cran()
	expect_output(greycurve(), NA)
	expect_output(greycurve(scalebar = TRUE), NA)
	expect_output(greycurve(1:100), NA)
})