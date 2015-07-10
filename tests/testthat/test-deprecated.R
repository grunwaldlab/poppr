context("Deprecated function tests")

data("microbov", package = "adegenet")

test_that("hierarchy methods all work, but return warnings", {
  skip_on_cran()
  # sethierarchy
  micstrat <- strata(microbov, value = data.frame(other(microbov)))
  expect_warning(oldstrat <- sethierarchy(microbov, data.frame(other(microbov))), 
                 "strata\\(.+?value")
  # replacement methods can't match the call
  expect_warning(sethierarchy(microbov) <- data.frame(other(microbov)),
                 "strata\\()")
  # Old methods behave in the same way.
  expect_that(oldstrat, equals(micstrat))
  expect_that(microbov, equals(micstrat))
  
  # namehierarchy
  micname <- nameStrata(microbov, ~Country/Breed/Species)
  expect_warning(oldname <- namehierarchy(microbov, ~Country/Breed/Species),
                 "nameStrata\\(microbov")
  expect_warning(namehierarchy(microbov) <- ~Country/Breed/Species,
                 "nameStrata\\()")
  expect_equal(oldname, micname)
  expect_equal(microbov, micname)
  
  # addhierarchy
  newstrat <- data.frame(TEST = sample(LETTERS, nInd(microbov), replace = TRUE))
  micadd <- addStrata(microbov, newstrat)
  expect_warning(oldadd <- addhierarchy(microbov, newstrat),
                 "addStrata\\(microbov")
  expect_warning(addhierarchy(microbov) <- newstrat,
                 "addStrata\\()")
  expect_equal(oldadd, micadd)
  expect_equal(microbov, micadd)

  # gethierarchy
  micget <- strata(microbov, ~Country/Breed/Species, combine = TRUE)
  micget2 <- strata(microbov, ~Country/Breed/Species, combine = FALSE)
  expect_warning(oldget <- gethierarchy(microbov, ~Country/Breed/Species),
                 "strata\\(microbov")
  expect_warning(oldget2 <- gethierarchy(microbov, ~Country/Breed/Species, combine = FALSE),
                 "strata\\(microbov.+?combine = FALSE")
  expect_equal(oldget, micget)
  expect_equal(oldget2, micget2)
  
  # splithierarchy
  strata(microbov) <- micget[3]
  micsplit <- splitStrata(microbov, ~Country/Breed/Species)
  expect_warning(oldsplit <- splithierarchy(microbov, ~Country/Breed/Species),
                 "splitStrata\\(microbov")
  expect_warning(splithierarchy(microbov) <- ~Country/Breed/Species,
                 "splitStrata\\()")
  expect_equal(oldsplit, micsplit)
  expect_equal(microbov, micsplit)
  expect_equal(microbov, micname)
})

test_that("setpop method works", {
  skip_on_cran()
  strata(microbov) <- data.frame(other(microbov))
  nameStrata(microbov) <- ~Country/Breed/Species

  micpop <- setPop(microbov, ~Country/Breed)
  expect_warning(oldpop <- setpop(microbov, ~Country/Breed),
                 "setPop\\(microbov")
  expect_warning(setpop(microbov) <- ~Country/Breed,
                 "setPop\\()")
  expect_equal(oldpop, micpop)
  expect_equal(microbov, micpop)
})