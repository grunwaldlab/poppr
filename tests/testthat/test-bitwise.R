context("Bitwise distance calculations")
set.seed(991)
dat <- sample(c(0, 1, NA), 50, replace = TRUE, prob = c(0.49, 0.5, 0.01))
mat <- matrix(dat, nrow = 5, ncol = 10)
mat2 <- mat
mat2[-1, ] <- mat2[-1, ] * 2

test_that("bitwise.dist works for haploids", {
  mat.gl <- new("genlight", mat, parallel = FALSE)
  expected <- rowSums(apply(mat, 2, dist), na.rm = TRUE)/10
  expect_equivalent(as.vector(bitwise.dist(mat.gl, threads = 1L)), expected)
  expect_true(all(ploidy(mat.gl) == 1))
})

test_that("bitwise.dist works for diploids with only one type of homozygote", {
  mat2.gl <- new("genlight", mat2, parallel = FALSE)
  expect_error(bitwise.dist(mat2.gl, threads = 1L), "ploidy")
  ploidy(mat2.gl) <- rep(2, 5)
  expected <- rowSums(apply(mat2, 2, dist), na.rm = TRUE)/20
  expect_equivalent(as.vector(bitwise.dist(mat2.gl, threads = 1L)), expected)
})

test_that("bitwise.dist can do euclidean", {
  mat2.gl <- new("genlight", mat2, parallel = FALSE)
  ploidy(mat2.gl) <- rep(2, 5)
  expect_equivalent(bitwise.dist(mat2.gl, scale_missing = TRUE, euclid = TRUE, threads = 1L), dist(mat2.gl))
})

test_that("bitwise.dist can do euclidean with lots of missing data", {
  # skip_on_cran()
  set.seed(999)
  mat2[sample(length(mat2), 10)] <- NA
  mat2.gl <- new("genlight", mat2, parallel = FALSE)
  ploidy(mat2.gl) <- rep(2, 5)
  expect_equivalent(bitwise.dist(mat2.gl, scale_missing = TRUE, euclid = TRUE, threads = 1L), dist(mat2.gl))
})

test_that("bitwise.dist can actually handle genind objects", {
  # skip_on_cran()
  data("partial_clone", package = "poppr")
  pdist <- diss.dist(partial_clone, percent = TRUE)
  expect_equivalent(pdist, bitwise.dist(partial_clone))
})


test_that("bitwise.dist produces reasonable results for diploids", {

  # skip_on_cran()
  dat <- list(c(2,2,2,2,2,2,2,2,2,0),
              c(1,1,1,0,0,0,0,0,0,2),
              c(2,2,2,2,2,2,2,2,2,2),
              c(2,2,2,2,2,2,2,2,2,0),
              c(2,NA,NA,NA,NA,NA,NA,NA,NA,NA))
  z <- new("genlight",dat, parallel = FALSE)
  
  # Sun May 13 12:37:19 2018 ------------------------------
  # Jonah originally set this up so it was simple to test, but did not account
  # for the fact that we needed to test when the missing data were in both the
  # i and the j positions. In this case, the missing are only in the j position.
  # Luckily, we can easily test this by inverting the data and results. In order
  # to distinguish between the two forms, I'll append i and j.
  missing_match_dif_i <- bitwise.dist(z[5:1, ], missing_match = TRUE, mat = TRUE, differences_only = TRUE)
  missing_match_dif_j <- bitwise.dist(z, missing_match = TRUE, mat = TRUE, differences_only = TRUE)
  
  expected_match_dif <- c(0.0, 1.0, 0.1, 0.0, 0.0, 
                          1.0, 0.0, 0.9, 1.0, 0.1, 
                          0.1, 0.9, 0.0, 0.1, 0.0, 
                          0.0, 1.0, 0.1, 0.0, 0.0, 
                          0.0, 0.1, 0.0, 0.0, 0.0)
  dim(expected_match_dif) <- c(5L, 5L)
  
  missing_nomatch_dif_i <- bitwise.dist(z[5:1, ], missing_match = FALSE, mat = TRUE, differences_only = TRUE)
  missing_nomatch_dif_j <- bitwise.dist(z, missing_match = FALSE, mat = TRUE, differences_only = TRUE)
  
  expected_nomatch_dif <- c(0.0, 1.0, 0.1, 0.0, 0.9, 
                            1.0, 0.0, 0.9, 1.0, 1.0, 
                            0.1, 0.9, 0.0, 0.1, 0.9, 
                            0.0, 1.0, 0.1, 0.0, 0.9, 
                            0.9, 1.0, 0.9, 0.9, 0.0)
  dim(expected_nomatch_dif) <- c(5L, 5L)
  
  expect_equivalent(missing_match_dif_i,   expected_match_dif[5:1, 5:1])
  expect_equivalent(missing_nomatch_dif_i, expected_nomatch_dif[5:1, 5:1])
  
  expect_equivalent(missing_match_dif_j,   expected_match_dif)
  expect_equivalent(missing_nomatch_dif_j, expected_nomatch_dif)
  
  missing_match_dist_i <- bitwise.dist(z[5:1, ], missing_match = TRUE, mat = TRUE, differences_only = FALSE)
  missing_match_dist_j <- bitwise.dist(z, missing_match = TRUE, mat = TRUE, differences_only = FALSE)
  expected_match_dist <- c(0.0, 17.0/20.0, 2.0/20.0, 0.0, 0.0, 
                           17.0/20.0, 0.0, 15.0/20.0, 17.0/20.0, 1.0/20.0, 
                           2.0/20.0, 15.0/20.0, 0.0, 2.0/20.0, 0.0, 
                           0.0, 17.0/20.0, 2.0/20.0, 0.0, 0.0, 
                           0.0, 1.0/20.0, 0.0, 0.0, 0.0)
  dim(expected_match_dist) <- c(5L, 5L)
  expect_equivalent(missing_match_dist_i, expected_match_dist[5:1, 5:1])
  expect_equivalent(missing_match_dist_j, expected_match_dist)
})

test_that("bitwise.ia produce reasonable results for haploids", {
  # skip_on_cran()
  dat <- list(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
              c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
              c(1, 1, NA, NA, NA, NA, NA, NA, NA, NA))
  z   <- new("genlight",dat, parallel = FALSE)
  expected_dist      <- c(0, 1, 0, 1, 0, 0.2, 0, 0.2, 0)
  dim(expected_dist) <- c(3L, 3L)
  
  missing_match_dif_i <- bitwise.dist(z[3:1, ], missing_match = TRUE, mat = TRUE, differences_only = TRUE)
  missing_match_dif_j <- bitwise.dist(z, missing_match = TRUE, mat = TRUE, differences_only = TRUE)
  
  expect_equivalent(missing_match_dif_i, expected_dist[3:1, 3:1])
  expect_equivalent(missing_match_dif_j, expected_dist)
  
  expected_nomatch_dist <- c(0, 1, .8, 1, 0, 1, .8, 1, 0)
  dim(expected_nomatch_dist) <- c(3L, 3L)

  missing_nomatch_dif_i <- bitwise.dist(z[3:1, ], missing_match = FALSE, mat = TRUE, differences_only = TRUE)
  missing_nomatch_dif_j <- bitwise.dist(z, missing_match = FALSE, mat = TRUE, differences_only = TRUE)
  
  expect_equivalent(missing_nomatch_dif_i, expected_nomatch_dist[3:1, 3:1])
  expect_equivalent(missing_nomatch_dif_j, expected_nomatch_dist)
  
})

context("bitwise.ia cromulence")

test_that("bitwise.ia can use both missing-match and missing-nomatch ", {
  # skip_on_cran()
  dat <- list(c(2,2,2,2,2,2,2,2,2,0),
              c(1,1,1,0,0,0,0,0,0,2),
              c(2,2,2,2,2,2,2,2,2,2),
              c(2,2,2,2,2,2,2,2,2,0),
              c(2,NA,NA,NA,NA,NA,NA,NA,NA,NA))
  z <- new("genlight",dat, parallel = FALSE)
  tz <- tab(z, NA.method = "asis")
  tz[tz == 2] <- "A/A"
  tz[tz == "1"] <- "A/B"
  tz[tz == "0"] <- "B/B"
  tzg <- df2genind(tz, sep = "/")
  
  # Missing-match is the default and matches that of the non-bitwise version
  expect_equal(bitwise.ia(z, missing_match = TRUE), ia(tzg)[[2]])
  expect_equal(bitwise.ia(z[5:1, ], missing_match = TRUE), ia(tzg)[[2]])
  
  ianm <- function(i, z){
    as.vector(bitwise.dist(z[, i], percent = FALSE, missing_match = FALSE))
  }
  np  <- choose(5, 2)
  bdl <- vapply(seq(10), ianm, integer(np), z)
  dlist <- list(
             d.vector = colSums(bdl),
             d2.vector = colSums(bdl * bdl),
             D.vector = rowSums(bdl)
           )
  res <- poppr:::ia_from_d_and_D(dlist, np)
  expect_equal(res[[2]], bitwise.ia(z, missing_match = FALSE))
})
