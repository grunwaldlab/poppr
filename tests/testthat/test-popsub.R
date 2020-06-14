context("Population subset tests")

test_that("subsetting needs a genind object", {
  expect_error(popsub(1:10), "popsub requires a genind or genlight object\n")
})

test_that("sublist needs to match populations", {
  data(microbov, package = "adegenet")
  expect_error(popsub(microbov, "missingno"), 
               poppr:::unmatched_pops_warning(popNames(microbov), "missingno"))
})

test_that("sum of the pops equal the whole", {
  data(microbov, package = "adegenet")
  mb1to7  <- popsub(microbov, 1:7)
  mb8to15 <- popsub(microbov, 8:15)
  expect_that(nInd(mb1to7) + nInd(mb8to15), equals(nInd(microbov)))
  })

test_that("subsetting works with populations", {
  data(nancycats, package = "adegenet")
  temp  <- nancycats@pop=="P04" | nancycats@pop=="P08"
  p48   <- nancycats[temp, ]
  p4    <- nancycats[nancycats@pop == "P04", , drop = TRUE]
  nan48 <- popsub(nancycats, c(4, 8), drop = FALSE)
  
  # Matrices equivalent
  expect_that(nan48@tab, equals(p48@tab))
  ## Dropping columns
  expect_that(popsub(nan48, "P04")@tab, equals(p4@tab))
  expect_that(popsub(nan48, 1)@tab, equals(p4@tab))
  expect_that(popsub(nancycats, 4)@tab, equals(p4@tab))
  # Populations equivalent
  expect_that(as.character(pop(nan48)), is_identical_to(as.character(pop(p48))))
  # Individuals equivalent
  expect_that(indNames(nan48), is_identical_to(indNames(p48)))
  # Rejects unknown populations
  expect_that(popsub(nancycats, 18), gives_warning())
  # Rejects old arguments
  expect_warning({
    popsub(nancycats, blacklist = 1:5)
  }, "exclude = 1:5")
  # Rejects equivalent exclude and sublist
  ## As numeric
  expect_that(popsub(nancycats, sublist = 1, exclude = 1), gives_warning())
  ## As characters
  expect_that(popsub(nancycats, sublist = "P01", exclude = "P01"), gives_warning())
  expect_that(popsub(nancycats, sublist = "P01", exclude = 1), gives_warning())
  expect_that(popsub(nancycats, sublist = 1, exclude = "P01"), gives_warning())
  expect_that(popsub(nancycats, sublist = c(4, 8), exclude = c(4, 8)), 
              gives_warning())
  expect_that(popsub(nancycats, sublist = c(4, 8), exclude = "P08")@tab, 
              equals(p4@tab))
  # numeric and character are the same
  expect_that(popsub(nancycats, c("P04", "P08"), drop = FALSE)@tab, equals(nan48@tab))
})

test_that("subsetting doesn't work without populations", {
  data(partial_clone, package = "poppr")
  p1 <- 1:50 %% 4 == 1
  expect_that(nInd(popsub(partial_clone, 1)), equals(nInd(partial_clone[p1, ])))
  expect_that(popsub(partial_clone, 1)@tab, equals(partial_clone[p1, , drop = TRUE]@tab))
  pop(partial_clone) <- NULL
  expect_that(popsub(partial_clone, 1), gives_warning())
})

test_that("subsetting works with genclone objects", {
  data(nancycats, package = "adegenet")
  temp      <- nancycats@pop=="P04" | nancycats@pop=="P08"
  p48       <- nancycats[temp, ]
  p4        <- nancycats[nancycats@pop == "P04", , drop = TRUE]
  nancycats <- as.genclone(nancycats)
  nan48     <- popsub(nancycats, c(4, 8), drop = FALSE)
  
  # Matrices equivalent
  expect_that(nan48@tab, equals(p48@tab))
  ## Dropping columns
  expect_that(popsub(nan48, "P04")@tab, equals(p4@tab))
  expect_that(popsub(nan48, 1)@tab, equals(p4@tab))
  expect_that(popsub(nancycats, 4)@tab, equals(p4@tab))
  # Populations equivalent
  expect_that(as.character(pop(nan48)), is_identical_to(as.character(pop(p48))))
  # Individuals equivalent
  expect_that(indNames(nan48), is_identical_to(indNames(p48)))
  # Rejects unknown populations
  expect_that(popsub(nancycats, 18), gives_warning())
  # Rejects equivalent exclude and sublist
  ## As numeric
  expect_that(popsub(nancycats, sublist = 1, exclude = 1), gives_warning())
  ## As characters
  expect_that(popsub(nancycats, sublist = "P01", exclude = "P01"), gives_warning())
  expect_that(popsub(nancycats, sublist = "P01", exclude = 1), gives_warning())
  expect_that(popsub(nancycats, sublist = 1, exclude = "P01"), gives_warning())
  expect_that(popsub(nancycats, sublist = c(4, 8), exclude = c(4, 8)), 
              gives_warning())
  expect_that(popsub(nancycats, sublist = c(4, 8), exclude = "P08")@tab, 
              equals(p4@tab))
  # numeric and character are the same
  expect_that(popsub(nancycats, c("P04", "P08"), drop = FALSE)@tab, equals(nan48@tab))
})
