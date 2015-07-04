context("Amova tests")

data("Aeut", package = "poppr")
strata(Aeut) <- other(Aeut)$population_hierarchy[-1]

agc    <- as.genclone(Aeut)
Athena <- popsub(agc, "Athena")
res    <- poppr.amova(Aeut, ~Pop/Subpop, quiet = TRUE)
rescc  <- poppr.amova(Aeut, ~Pop/Subpop, quiet = TRUE, clonecorrect = TRUE)


test_that("Amova returns published values", {
	expect_that(res$componentsofcovariance[, 2], equals(c(70.0067859292295, 
														  8.40748251295027, 
														  21.5857315578203, 
														  100)))
	expect_that(res$componentsofcovariance[, 1], equals(c(11.0634458464745, 
														  1.3286673034988, 
														  3.41127747798475, 
														  15.8033906279581)))
	expect_that(rescc$componentsofcovariance[, 2], equals(c(66.7776803325885, 
															6.10535452127678, 
															27.1169651461347, 
															100)))
	expect_that(rescc$componentsofcovariance[, 1], equals(c(10.4131525344064, 
															0.952054452775842, 
															4.22855500416477, 
															15.593761991347)))
	})
test_that("AMOVA handles subsetted genclone objects", {
	Athena.mlg <- mlg.vector(Athena)
	agc.mlg <- mlg.vector(agc)
	Athena.AMOVA <- poppr.amova(Athena, ~Subpop, quiet = TRUE)

	# All MLGs are represented
	expect_false(anyNA(match(1:max(agc.mlg), agc.mlg)))

	# Some MLGS are missing
	expect_true(anyNA(match(1:max(Athena.mlg), Athena.mlg)))

	#AMOVA will still work on the data set with missing MLGs
	expect_that(class(Athena.AMOVA), equals("amova"))	
	})


test_that("AMOVA can take a distance matrix as input", {
  skip_on_cran()
  adist   <- diss.dist(clonecorrect(Aeut, strata = NA), percent = FALSE)
  root_adist  <- sqrt(adist)
  
  # UNSQUARED
  usqres <- poppr.amova(Aeut, ~Pop/Subpop, quiet = TRUE, dist = root_adist, 
                         squared = FALSE)
  # SQUARED
  sqres <- poppr.amova(Aeut, ~Pop/Subpop, quiet = TRUE, dist = adist, 
                        squared = TRUE)
  
  expect_equivalent(usqres$componentsofcovariance, res$componentsofcovariance)
  expect_equivalent(usqres$componentsofcovariance, sqres$componentsofcovariance)
  err <- c("Uncorrected.+?97.+?Clone.+?70.+?provided.+?119")
  expect_error(poppr.amova(Athena, ~Subpop, quiet = TRUE, dist = adist), err)
})


test_that("AMOVA can work on clone correction from mlg.filter and do filtering", {
  skip_on_cran()
  data("monpop", package = "poppr")
  splitStrata(monpop) <- ~Tree/Year/Symptom
  expect_warning(poppr.amova(monpop, ~Symptom/Year), "Zero")
  mondist <- diss.dist(monpop)
  monfilt <- mlg.filter(monpop, distance = mondist, threshold = 1.1, 
                        stats = "DIST")
  mlg.filter(monpop, distance = mondist) <- 1.1 
  monfiltdist <- as.dist(monfilt)
  monres <- poppr.amova(monpop, ~Symptom/Year, dist = monfiltdist)
  msg <- "Original.+?264"
  expect_message(res <- poppr.amova(monpop, ~Symptom/Year, filter = TRUE, threshold = 1.1), msg)
  expect_equivalent(monres$componentsofcovariance, res$componentsofcovariance)
})

test_that("AMOVA can calculate within individual variance for diploids", {
  skip_on_cran()
  data("microbov", package = "adegenet")
  strata(microbov) <- data.frame(other(microbov))
  mics <- microbov[pop = 1:2]
  expect_output(micwithin  <- poppr.amova(mics, ~breed), "Removing")
  expect_output(poppr.amova(mic, ~breed, correction = "cai"), "Cailliez constant")
  expect_output(poppr.amova(mic, ~breed, correction = "lin"), "Lingoes constant")
  expect_output(micwithout <- poppr.amova(mics, ~breed, within = FALSE), "Removing")
  expect_equal(dim(micwithin$componentsofcovariance), c(4, 2))
  expect_equal(dim(micwithout$componentsofcovariance), c(3, 2))
})