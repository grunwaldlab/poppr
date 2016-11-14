context("bootclass methods")

data("nancycats", package = "adegenet")
nan9 <- nancycats[pop = 9]
nanreps <- structure(c(1.99999, 2, 1.99999, 2, 2, 2, 1.99999, 1.99999, 2
), .Names = c("fca8", "fca23", "fca43", "fca45", "fca77", "fca78", 
"fca90", "fca96", "fca37"))
bgnan <- new("bootgen", nan9, na = "asis", freq = FALSE)
bvnan <- new("bruvomat", nan9, nanreps)

test_that("bootgen can be subset correctly without loci specified", {
  expect_identical(bgnan[], bgnan)
})

test_that("bruvomat can be subset correctly without loci specified", {
  expect_equivalent(bvnan[], bvnan)
})

test_that("bruvomat can be subset correctly without loci specified", {
  subafter <- bvnan@mat[, rep((1:9) %% 2 == 1, each = 2)]
  subfore  <- bvnan[, (1:9) %% 2 == 1]@mat
  dimnames(subafter) <- NULL -> dimnames(subfore)
  expect_equivalent(subfore, subafter)
})