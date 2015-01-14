context("Genclone coercion tests")

test_that("A genclone object contains a genind object", {

	data(partial_clone, package = "poppr")
	pc <- as.genclone(partial_clone)
	expect_that(slotNames(partial_clone), equals(slotNames(pc)[-c(1:2)]))
	expect_that(partial_clone@tab, equals(pc@tab))
	expect_that(partial_clone@loc.names, is_identical_to(pc@loc.names))
	expect_that(partial_clone@loc.nall, is_identical_to(pc@loc.nall))
	expect_that(partial_clone@all.names, is_identical_to(pc@all.names))
	expect_that(partial_clone@ind.names, is_identical_to(pc@ind.names))
	expect_that(partial_clone@pop, is_identical_to(pc@pop))
	expect_that(partial_clone@pop.names, is_identical_to(pc@pop.names))
	expect_that(partial_clone@ploidy, is_identical_to(pc@ploidy))
	expect_that(partial_clone@type, is_identical_to(pc@type))
	expect_that(partial_clone@other, is_identical_to(pc@other))
	expect_that(pc@mlg[], is_identical_to(mlg.vector(partial_clone)))
})

test_that("Hierarchy methods work for genclone objects.", {
  data(Aeut, package = "poppr")
  agc <- as.genclone(Aeut)
  expect_that(length(gethierarchy(agc)), equals(3))
  expect_that(agc@pop.names, equals(c(P1 = "Athena", P2 = "Mt. Vernon")))
  expect_that({agcsplit <- splithierarchy(agc, ~Pop/Subpop)}, gives_warning())
  expect_that(gethierarchy(agcsplit), equals(gethierarchy(agc, ~Pop/Subpop, combine = FALSE)))
  expect_that(sethierarchy(agc, gethierarchy(agcsplit)), equals(agcsplit))
  namehierarchy(agcsplit) <- ~Field/Core
  expect_that(names(gethierarchy(agcsplit)), equals(c("Field", "Core")))
  setpop(agc) <- ~Pop/Subpop
  expect_that(agc@pop.names, equals(c("Athena_1", "Athena_2", "Athena_3", 
                                      "Athena_4", "Athena_5", "Athena_6", 
                                      "Athena_7", "Athena_8", "Athena_9", 
                                      "Athena_10", "Mt. Vernon_1", 
                                      "Mt. Vernon_2", "Mt. Vernon_3", 
                                      "Mt. Vernon_4", "Mt. Vernon_5", 
                                      "Mt. Vernon_6", "Mt. Vernon_7", 
                                      "Mt. Vernon_8")))
})