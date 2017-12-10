context("Genclone coercion tests")

data("Aeut", package = "poppr")
agc         <- as.genclone(Aeut)
strata(agc) <- other(agc)$population_hierarchy

data(partial_clone, package = "poppr")
pc <- as.genclone(partial_clone)

data("old_Pinf", package = "poppr")
data("old_partial_clone", package = "poppr")

data("Pinf", package = "poppr")

test_that("A genclone object contains a genind object", {
	expect_equal(slotNames(partial_clone), slotNames(pc)[-1])
	expect_equal(tab(partial_clone), tab(pc))
	expect_identical(nAll(partial_clone), nAll(pc))
	expect_equivalent(alleles(partial_clone), alleles(pc))
	expect_identical(pop(partial_clone), pop(pc))
	expect_identical(popNames(partial_clone), popNames(pc))
	expect_identical(ploidy(partial_clone), ploidy(pc))
	expect_identical(partial_clone@type, pc@type)
	expect_identical(partial_clone@other, pc@other)
	expect_identical(pc@mlg[], mlg.vector(partial_clone))
})

test_that("A genclone object is a valid object", {
  skip_on_cran()
  expect_true(validObject(pc))
  expect_true(poppr:::valid.genclone(pc))
  pcno <- pc
  pcno@mlg <- 1:26
  expect_error(validObject(pcno), gettext("invalid class", domain = "R"))
  expect_message(poppr:::valid.genclone(pcno), 
                 gettext("Multilocus genotypes do not match the number of observations", domain = "R"))
})

test_that("Strata methods work for genclone objects.", {
  
  expect_equal(length(strata(agc)), 3)
  expect_equal(popNames(agc), c("Athena", "Mt. Vernon"))
  expect_that({agcsplit <- splitStrata(agc, ~Pop/Subpop)}, gives_warning())
  expect_equal(strata(agcsplit), strata(agc, ~Pop/Subpop, combine = FALSE))
  expect_equal(strata(agc, value = strata(agcsplit)), agcsplit)
  nameStrata(agcsplit) <- ~Field/Core
  expect_equal(names(strata(agcsplit)), c("Field", "Core"))
  setPop(agc) <- ~Pop/Subpop
  expect_that(popNames(agc), equals(c("Athena_1", "Athena_2", "Athena_3", 
                                      "Athena_4", "Athena_5", "Athena_6", 
                                      "Athena_7", "Athena_8", "Athena_9", 
                                      "Athena_10", "Mt. Vernon_1", 
                                      "Mt. Vernon_2", "Mt. Vernon_3", 
                                      "Mt. Vernon_4", "Mt. Vernon_5", 
                                      "Mt. Vernon_6", "Mt. Vernon_7", 
                                      "Mt. Vernon_8")))
})

test_that("Subsetting works with populations", {
  sA <- Pinf[pop = "South America"]
  nA <- Pinf[pop = "North America"]
  msg <- '\\) <- "original"'
  expect_warning(msA <- mll(sA), paste0("mll\\(sA", msg))
  expect_warning(mnA <- mll(nA), paste0("mll\\(nA", msg))
  expect_warning(mP  <- mll(Pinf), paste0("mll\\(Pinf", msg))
  expect_true(all(sort(c(msA, mnA)) == sort(mP)))
})

test_that("genclone2genind conversion works", {
  AGC <- genclone2genind(agc)
  expect_false(is.clone(AGC))
  expect_true(is.clone(as.genclone(AGC)))
})

test_that("old2new_genclone conversion works",{
  expect_is(old2new_genclone(old_Pinf), "genclone")
  expect_is(old2new_genclone(old_partial_clone), "genind")
})

test_that("print method will show all populations", {
  setPop(agc) <- ~Pop/Subpop
  expect_output(print(agc), paste(popNames(agc)[1:5], collapse = " "))
  expect_output(print(agc, fullnames = FALSE), "...")
})

test_that("subsetting a genclone object retains the MLG definitions", {
  idn <- mlg.id(pc)[[1]]
  idl <- indNames(pc) %in% idn
  idi <- which(idl)
  expect_equal(mll(pc[idn]), mll(pc[idl]))
  expect_equal(mll(pc[idi]), mll(pc[idl]))
})

test_that("genclone objects can be subset with character strings or factors", {
  idn <- mlg.id(pc)[[1]]
  idf <- as.factor(idn)
  expect_equal(mll(pc[idn]), mll(pc[idf]))
})

test_that("mlg.filter<- works for genclone, but not genind", {
  skip_on_cran()
  mlg.filter(pc, distance = nei.dist) <- 0.1
  expect_equal(nmll(pc), 24)
  expect_warning(mlg.filter(partial_clone) <- 0.1)
})