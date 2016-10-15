context("Test trees")

library("ape")
data("nancycats", package = "adegenet")
data("Aeut", package = "poppr")

myFastME     <- function(x) fastme.bal(x, nni = TRUE, spr = FALSE, tbr = TRUE)

nan9         <- popsub(nancycats, 9)
nanreps      <- fix_replen(nancycats, rep(2, 9))

strata(Aeut) <- other(Aeut)$population_hierarchy[-1]
Aeut.pop     <- genind2genpop(Aeut, pop = ~Pop/Subpop, quiet = TRUE)

test_that("bruvo.boot produces a tree from any tree function", {
	skip_on_cran()
	nanupgma <- bruvo.boot(nan9, replen = nanreps, sample = 20, quiet = TRUE)
	nannj    <- bruvo.boot(nan9, replen = nanreps, sample = 20, quiet = TRUE, tree = "nj")
	nanfast  <- bruvo.boot(nan9, replen = nanreps, sample = 20, quiet = TRUE, tree = myFastME)

	expect_true(ape::is.ultrametric(nanupgma))
	expect_false(ape::is.ultrametric(nannj))
	expect_false(ape::is.ultrametric(nanfast))
})

test_that("bruvo.boot rejects non-ssr data", {
	expect_error(bruvo.boot(Aeut))
})

test_that("aboot can take in matrices", {
  skip_on_cran()
  Aeut.tab <- tab(Aeut.pop)
  apop <- aboot(Aeut.pop, showtree = FALSE, quiet = TRUE)
  amat <- aboot(Aeut.tab, showtree = FALSE, quiet = TRUE)
  expect_equal(apop$tip.label, amat$tip.label)
  dimnames(Aeut.tab) <- NULL
  amat <- aboot(Aeut.tab, quiet = TRUE)
  expect_equal(amat$tip.label, .genlab("", nrow(Aeut.tab)))
})

test_that("aboot can also produce trees from any function", {
	skip_on_cran()

	nanupgma <- aboot(nan9, dist = rogers.dist, sample = 20, quiet = TRUE)
	nannj    <- aboot(nan9, dist = rogers.dist, sample = 20, quiet = TRUE, tree = "nj")
	nanfast  <- aboot(nan9, dist = rogers.dist, sample = 20, quiet = TRUE, tree = myFastME)

	expect_true(ape::is.ultrametric(nanupgma))
	expect_false(ape::is.ultrametric(nannj))
	expect_false(ape::is.ultrametric(nanfast))
})

test_that("aboot can utilize anonymous functions", {
	skip_on_cran()
	nantree <- aboot(nan9, dist = function(x) dist(x$tab, method = "manhattan"), sample = 20, quiet = TRUE)
	expect_is(nantree, "phylo")
})

test_that("aboot can utilize bootgen2genind", {
  skip_on_cran()
  suppressWarnings(nantreebg <- aboot(nan9, dist = function(x) nei.dist(bootgen2genind(x)), sample = 20, quiet = TRUE))
  suppressWarnings(nantree   <- aboot(nan9, dist = "nei.dist", sample = 20, quiet = TRUE))
  expect_is(nantreebg, "phylo")
  expect_is(nantree, "phylo")
  expect_equivalent(nantreebg$edge.length, nantree$edge.length)
})


test_that("aboot can handle populations", {
  skip_on_cran()
	set.seed(999)
	expect_warning(Atreen <- aboot(Aeut.pop, tree = "nj", sample = 20, quiet = TRUE), "Felsenstein")
	Atreeu                <- aboot(Aeut.pop, sample = 20, quiet = TRUE)
	expect_warning(AtreeF <- aboot(Aeut.pop, tree = myFastME, sample = 20, quiet = TRUE), "Felsenstein")

	expect_true(ape::is.ultrametric(Atreeu))
	expect_false(ape::is.ultrametric(Atreen))
	expect_false(ape::is.ultrametric(AtreeF))
})

test_that("aboot can convert genind to genpop", {
  skip_on_cran()
  set.seed(999)
  Atree <- aboot(Aeut, strata = ~Pop/Subpop, sample = 20, quiet = TRUE)
  
  expect_equal(ape::Ntip(Atree), 18L)
})

test_that("aboot can handle genlight objects", {
	skip_on_cran()
	set.seed(999)
	gc <- as.snpclone(glSim(100, 0, n.snp.struc = 1e3, ploidy = 2, parallel = FALSE))
	gtree <- aboot(gc, distance = bitwise.dist, sample = 5, quiet = TRUE)
	expect_warning(gwarn <- aboot(gc, strata = ~none, distance = bitwise.dist, sample = 5, quiet = TRUE))
	expect_is(gtree, "phylo")
	expect_is(gwarn, "phylo")
	expect_equal(ape::Ntip(gtree), ape::Ntip(gwarn))
})
