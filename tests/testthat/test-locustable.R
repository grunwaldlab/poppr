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
	randall              <- sample(nLoc(nancy), 1)
	expect_message(nanlt <- locus_table(nancy), "Simpson")
	expect_message(nanlt <- locus_table(nancy, index = "shannon"), "Shannon")
	expect_message(nanlt <- locus_table(nancy, index = "invsimpson"), "Taylor")
	expect_message(nangt <- locus_table(nancy, lev = "genotype"), "genotype")
	expect_gt(nangt[randall, "genotype"], nanlt[randall, "allele"])
	expect_output(nanlt <- locus_table(nancy, information = FALSE), NA)
})

test_that("locus_table will accurately calculate Hexp", {
  skip_on_cran()
  
  # From Kosman, 2003: http://onlinelibrary.wiley.com/doi/10.1046/j.1365-3059.2003.00923.x/full
  qs <- c(2/3, 1, 5/6, 2/3, 1/2, 1/3, 1/6, 1/6)
  hs <- c(4/9, 0, 5/18, 4/9, 1/2, 4/9, 5/18, 5/18)
  
  # Original table, haploid.
  xh <- "
  1 1 1 1 0 0 0 0
  0 1 1 1 1 1 1 0
  1 1 1 1 1 0 0 0
  0 1 1 1 1 1 0 0
  1 1 1 0 0 0 0 0
  1 1 0 0 0 0 0 1"
  xh_tab <- read.table(text = xh, sep = "")
  
  # Double the data, but same result for Hs because of allele frequencies
  xd_tab <- apply(xh_tab, 2, function(x) paste(x + 1, sample(x + 1), sep = "/"))
  xhf <- file()
  xdf <- file()
  write.table(xh_tab + 1, file = xhf, sep = ",")
  write.table(xd_tab, file = xdf, sep = ",")
  x_hap <- pegas::read.loci(xhf, loci.sep = ",")
  x_dip <- pegas::read.loci(xdf, loci.sep = ",", allele.sep = "/")
  close(xhf)
  close(xdf)
  
  lt_hap <- poppr::locus_table(pegas::loci2genind(x_hap))
  lt_dip <- poppr::locus_table(pegas::loci2genind(x_dip))
  
  expect_equivalent(lt_hap[1:8, "1-D"], hs)
  expect_equivalent(lt_dip[1:8, "1-D"], hs)

  expect_equivalent(lt_hap["mean", "1-D"], mean(hs))
  expect_equivalent(lt_dip["mean", "1-D"], mean(hs))
  
  expect_equivalent(lt_hap[1:8, "Hexp"], (6/5)*hs)
  expect_equivalent(lt_dip[1:8, "Hexp"], (12/11)*hs)
  
})