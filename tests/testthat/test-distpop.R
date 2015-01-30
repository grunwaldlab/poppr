context("Population Distance Tests")

test_that("dist.genpop matches distance", {
	data(nancycats, package = "adegenet")
	nanpop <- genind2genpop(nancycats, quiet = TRUE)
	nei      <- as.vector(nei.dist(nanpop))
	edwards  <- as.vector(edwards.dist(nanpop))
	reynolds <- as.vector(reynolds.dist(nanpop))
	rogers   <- as.vector(rogers.dist(nanpop))
	provesti <- as.vector(provesti.dist(nanpop))
	expect_that(as.vector(dist.genpop(nanpop, method = 1)), equals(nei))
	expect_that(as.vector(dist.genpop(nanpop, method = 2)), equals(edwards))
	expect_that(as.vector(dist.genpop(nanpop, method = 3)), equals(reynolds))
	expect_that(as.vector(dist.genpop(nanpop, method = 4)), equals(rogers))
	expect_that(as.vector(dist.genpop(nanpop, method = 5)), equals(provesti))
})

test_that("aboot works with diss.dist", {
  data(nancycats, package = "adegenet")
  nan1 <- popsub(nancycats, 9)
  ab <- aboot(nan1, dist = diss.dist, sample = 2, quiet = TRUE, showtree = FALSE)
  expect_that(class(ab), equals("phylo"))
})