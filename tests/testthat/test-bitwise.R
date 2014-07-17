context("Bitwise distance calculations")

test_that("bitwise_distance.c produces reasonable results", {
  dat <- list(c(2,2,2,2,2,2,2,2,2,0),c(1,1,1,0,0,0,0,0,0,2),c(2,2,2,2,2,2,2,2,2,2),c(2,2,2,2,2,2,2,2,2,0),c(2,NA,NA,NA,NA,NA,NA,NA,NA,NA))
  z <- new("genlight",dat)
  missing_match <- bitwise.dist(z,missing_match=TRUE,mat=TRUE)
  dim(missing_match) <- NULL
  expected_match <- c(0.0, 1.0, 0.1, 0.0, 0.0, 1.0, 0.0, 0.9, 1.0, 0.1, 0.1, 0.9, 0.0, 0.1, 0.0, 0.0, 1.0, 0.1, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0)
  missing_nomatch <- bitwise.dist(z,missing_match=FALSE,mat=TRUE)
  dim(missing_nomatch) <- NULL
  expected_nomatch <- c(0.0, 1.0, 0.1, 0.0, 0.9, 1.0, 0.0, 0.9, 1.0, 1.0, 0.1, 0.9, 0.0, 0.1, 0.9, 0.0, 1.0, 0.1, 0.0, 0.9, 0.9, 1.0, 0.9, 0.9, 0.0)
  expect_that(missing_match, is_equivalent_to(expected_match))
  expect_that(missing_nomatch, is_equivalent_to(expected_nomatch))
 })
