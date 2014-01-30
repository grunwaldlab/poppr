context("Population subset")

test_that("subsetting needs a genind object", {
  expect_error(popsub(1:10), "popsub requires a genind object\n")
})

test_that("sublist needs to match populations", {
  data(microbov)
  expect_error(popsub(microbov, "missingno"), 
               poppr:::unmatched_pops_warning(microbov@pop.names, "missingno"))
})

test_that("sum of the pops equal the whole", {
  data(microbov)
  mb1to7  <- popsub(microbov, 1:7)
  mb8to15 <- popsub(microbov, 8:15)
  expect_that(nInd(mb1to7) + nInd(mb8to15), equals(nInd(microbov)))
  })

test_that("subsetting works with populations", {
  data(nancycats)
  temp  <- nancycats@pop=="P04" | nancycats@pop=="P08"
  p48   <- nancycats[temp, ]
  p4    <- nancycats[nancycats@pop == "P04", , drop = TRUE]
  nan48 <- popsub(nancycats, c(4, 8), drop = FALSE)
  
  # Matrices equivalent
  expect_that(nan48@tab, equals(p48@tab))
  ## Dropping columns
  expect_that(popsub(nan48, "4")@tab, equals(p4@tab))
  expect_that(popsub(nan48, 1)@tab, equals(p4@tab))
  expect_that(popsub(nancycats, 4)@tab, equals(p4@tab))
  # Populations equivalent
  expect_that(as.character(pop(nan48)), matches(as.character(pop(p48))))
  # Individuals equivalent
  expect_that(nan48@ind.names, matches(p48@ind.names))
  # Rejects unknown populations
  expect_that(popsub(nancycats, 18), gives_warning())
  # Rejects equivalent blacklist and sublist
  ## As numeric
  expect_that(popsub(nancycats, sublist = 1, blacklist = 1), gives_warning())
  ## As characters
  expect_that(popsub(nancycats, sublist = "1", blacklist = "1"), gives_warning())
  expect_that(popsub(nancycats, sublist = "1", blacklist = 1), gives_warning())
  expect_that(popsub(nancycats, sublist = 1, blacklist = "1"), gives_warning())
  expect_that(popsub(nancycats, sublist = c(4, 8), blacklist = c(4, 8)), 
              gives_warning())
  expect_that(popsub(nancycats, sublist = c(4, 8), blacklist = "8")@tab, 
              equals(p4@tab))
  # numeric and character are the same
  expect_that(popsub(nancycats, c("4", "8"), drop = FALSE)@tab, equals(nan48@tab))
})

test_that("subsetting doesn't work without populations", {
  data(partial_clone)
  p1 <- 1:50 %% 4 == 1
  expect_that(nInd(popsub(partial_clone, 1)), equals(nInd(partial_clone[p1, ])))
  expect_that(popsub(partial_clone, 1)@tab, equals(partial_clone[p1, , drop = TRUE]@tab))
  pop(partial_clone) <- NULL
  expect_that(popsub(partial_clone, 1), gives_warning())
})

test_that("subsetting works with genclone objects", {
  data(nancycats)
  temp      <- nancycats@pop=="P04" | nancycats@pop=="P08"
  p48       <- nancycats[temp, ]
  p4        <- nancycats[nancycats@pop == "P04", , drop = TRUE]
  nancycats <- as.genclone(nancycats)
  nan48     <- popsub(nancycats, c(4, 8), drop = FALSE)
  
  # Matrices equivalent
  expect_that(nan48@tab, equals(p48@tab))
  ## Dropping columns
  expect_that(popsub(nan48, "4")@tab, equals(p4@tab))
  expect_that(popsub(nan48, 1)@tab, equals(p4@tab))
  expect_that(popsub(nancycats, 4)@tab, equals(p4@tab))
  # Populations equivalent
  expect_that(as.character(pop(nan48)), matches(as.character(pop(p48))))
  # Individuals equivalent
  expect_that(nan48@ind.names, matches(p48@ind.names))
  # Rejects unknown populations
  expect_that(popsub(nancycats, 18), gives_warning())
  # Rejects equivalent blacklist and sublist
  ## As numeric
  expect_that(popsub(nancycats, sublist = 1, blacklist = 1), gives_warning())
  ## As characters
  expect_that(popsub(nancycats, sublist = "1", blacklist = "1"), gives_warning())
  expect_that(popsub(nancycats, sublist = "1", blacklist = 1), gives_warning())
  expect_that(popsub(nancycats, sublist = 1, blacklist = "1"), gives_warning())
  expect_that(popsub(nancycats, sublist = c(4, 8), blacklist = c(4, 8)), 
              gives_warning())
  expect_that(popsub(nancycats, sublist = c(4, 8), blacklist = "8")@tab, 
              equals(p4@tab))
  # numeric and character are the same
  expect_that(popsub(nancycats, c("4", "8"), drop = FALSE)@tab, equals(nan48@tab))
})
