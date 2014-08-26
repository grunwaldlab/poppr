context("Bitwise distance and pgen calculations")

test_that("bitwise.dist produces reasonable results", {

# Required to circumvent a windows specific error in adegenet
  if ((.Platform)$OS.type == "windows"){
    mclapply <- function (X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE,
      mc.silent = FALSE, mc.cores = 1L, mc.cleanup = TRUE,
      mc.allow.recursive = TRUE){
        cores <- as.integer(mc.cores)
        if (cores < 1L)
          stop("'mc.cores' must be >= 1")
        if (cores > 1L)
          lapply(X, FUN, ...)
      }
  }

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

test_that("pgen produces reasonable results", {

# Required to circumvent a windows specific error in adegenet
  if ((.Platform)$OS.type == "windows"){
    mclapply <- function (X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE,
      mc.silent = FALSE, mc.cores = 1L, mc.cleanup = TRUE,
      mc.allow.recursive = TRUE){
        cores <- as.integer(mc.cores)
        if (cores < 1L)
          stop("'mc.cores' must be >= 1")
        if (cores > 1L)
          lapply(X, FUN, ...)
      }
  }

  single_pop <- new("genlight", list(0, 1, 1, 2), ploidy=c(2,2,2,2))
  results_single_pop <- pgen(single_pop,log=FALSE)
  expected_single_pop <- structure(c(0.25, 0.5, 0.5, 0.25), .Dim = c(4L, 1L))

  multi_pop <- new("genlight", list(0, 1, 1, 2), ploidy=c(2,2,2,2))
  multi_pop$pop <- as.factor(1:4)
  results_multi_pop <- pgen(multi_pop,log=FALSE)
  expected_multi_pop <- structure(c(1, 0.5, 0.5, 1), .Dim = c(4L, 1L))
 
  expect_that(results_single_pop, is_equivalent_to(expected_single_pop))
  expect_that(results_multi_pop, is_equivalent_to(expected_multi_pop))
  expect_that(pgen(multi_pop,log=FALSE,by.pop=FALSE), is_equivalent_to(expected_single_pop))
 })
