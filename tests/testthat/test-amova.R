context("Amova tests")

test_that("Amova returns published values", {
	data(Aeut, package = "poppr")
	res   <- poppr.amova(Aeut, ~Pop/Subpop, quiet = TRUE)
	rescc <- poppr.amova(Aeut, ~Pop/Subpop, quiet = TRUE, clonecorrect = TRUE)
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
	data(Aeut, package = "poppr")
	agc <- as.genclone(Aeut)
	Athena <- popsub(agc, "Athena")
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
