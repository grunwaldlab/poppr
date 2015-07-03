context("Population Distance Tests")

test_that("dist.genpop matches distance", {
	data(nancycats, package = "adegenet")
	nanpop   <- genind2genpop(nancycats, quiet = TRUE)
	nei      <- as.vector(nei.dist(nanpop))
	edwards  <- as.vector(edwards.dist(nanpop))
	reynolds <- as.vector(reynolds.dist(nanpop))
	rogers   <- as.vector(rogers.dist(nanpop))
	provesti <- as.vector(provesti.dist(nanpop))
	expect_equivalent(as.vector(dist.genpop(nanpop, method = 1)), nei)
	expect_equivalent(as.vector(dist.genpop(nanpop, method = 2)), edwards)
	expect_equivalent(as.vector(dist.genpop(nanpop, method = 3)), reynolds)
	expect_equivalent(as.vector(dist.genpop(nanpop, method = 4)), rogers)
	expect_equivalent(as.vector(dist.genpop(nanpop, method = 5)), provesti)
})

test_that("aboot works with diss.dist", {
  data(nancycats, package = "adegenet")
  nan1 <- popsub(nancycats, 9)
  ab <- aboot(nan1, dist = diss.dist, sample = 2, quiet = TRUE, showtree = FALSE)
  expect_equal(class(ab), "phylo")
})