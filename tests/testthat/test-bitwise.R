context("Bitwise distance and pgen calculations")
set.seed(991)
dat <- sample(c(0, 1, NA), 50, replace = TRUE, prob = c(0.49, 0.5, 0.01))
mat <- matrix(dat, nrow = 5, ncol = 10)

test_that("bitwise.dist works for haploids", {
  mat.gl <- new("genlight", mat, parallel = FALSE)
  expected <- rowSums(apply(mat, 2, dist), na.rm = TRUE)/10
  expect_equivalent(as.vector(bitwise.dist(mat.gl)), expected)
  expect_true(all(ploidy(mat.gl) == 1))
})

test_that("bitwise.dist works for diploids with only one type of homozygote", {
  mat2 <- mat
  mat2[-1, ] <- mat2[-1, ] * 2
  mat2.gl <- new("genlight", mat2, parallel = FALSE)
  expect_error(bitwise.dist(mat2.gl), "ploidy")
  ploidy(mat2.gl) <- rep(2, 5)
  expected <- rowSums(apply(mat2, 2, dist), na.rm = TRUE)/20
  expect_equivalent(as.vector(bitwise.dist(mat2.gl)), expected)
})


test_that("bitwise.dist can actually handle genind objects", {
  skip_on_cran()
  data("partial_clone", package = "poppr")
  pdist <- diss.dist(partial_clone, percent = TRUE)
  expect_equivalent(pdist, bitwise.dist(partial_clone))
})


test_that("bitwise.dist produces reasonable results", {

  skip_on_cran()
  dat <- list(c(2,2,2,2,2,2,2,2,2,0),
              c(1,1,1,0,0,0,0,0,0,2),
              c(2,2,2,2,2,2,2,2,2,2),
              c(2,2,2,2,2,2,2,2,2,0),
              c(2,NA,NA,NA,NA,NA,NA,NA,NA,NA))
  z <- new("genlight",dat, parallel = FALSE)
  
  missing_match_dif <- bitwise.dist(z, missing_match=TRUE, mat=TRUE, differences_only=TRUE)
  dim(missing_match_dif) <- NULL
  expected_match_dif <- c(0.0, 1.0, 0.1, 0.0, 0.0, 1.0, 0.0, 0.9, 1.0, 0.1, 0.1, 0.9, 0.0, 0.1, 0.0, 0.0, 1.0, 0.1, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0)
  missing_nomatch_dif <- bitwise.dist(z,missing_match=FALSE,mat=TRUE, differences_only=TRUE)
  dim(missing_nomatch_dif) <- NULL
  expected_nomatch_dif <- c(0.0, 1.0, 0.1, 0.0, 0.9, 1.0, 0.0, 0.9, 1.0, 1.0, 0.1, 0.9, 0.0, 0.1, 0.9, 0.0, 1.0, 0.1, 0.0, 0.9, 0.9, 1.0, 0.9, 0.9, 0.0)
  
  expect_equivalent(missing_match_dif, expected_match_dif)
  expect_equivalent(missing_nomatch_dif, expected_nomatch_dif)
  
  missing_match_dist <- bitwise.dist(z,missing_match=TRUE,mat=TRUE,differences_only=FALSE)
  dim(missing_match_dist) <- NULL
  expected_match_dist <- c(0.0, 17.0/20.0, 2.0/20.0, 0.0, 0.0, 17.0/20.0, 0.0, 15.0/20.0, 17.0/20.0, 1.0/20.0, 2.0/20.0, 15.0/20.0, 0.0, 2.0/20.0, 0.0, 0.0, 17.0/20.0, 2.0/20.0, 0.0, 0.0, 0.0, 1.0/20.0, 0.0, 0.0, 0.0)
  #missing_nomatch <- bitwise.dist(z,missing_match=FALSE,mat=TRUE,differences_only=FALSE)
  #dim(missing_nomatch) <- NULL
  #expected_nomatch <- c(0.0, 1.0, 0.1, 0.0, 0.9, 1.0, 0.0, 0.9, 1.0, 1.0, 0.1, 0.9, 0.0, 0.1, 0.9, 0.0, 1.0, 0.1, 0.0, 0.9, 0.9, 1.0, 0.9, 0.9, 0.0)
  
  expect_equivalent(missing_match_dist, expected_match_dist)
  #expect_equivalent(missing_nomatch, expected_nomatch)
 })
