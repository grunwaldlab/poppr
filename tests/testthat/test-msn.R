context("Minimum spanning network tests")

data("partial_clone")
pc <- as.genclone(partial_clone)
mll.custom(pc) <- LETTERS[mll(pc)]
mll.levels(pc)[mll.levels(pc) == "Q"] <- "M"
mll(pc) <- "custom"

test_that("msn works with custom MLLs", {
  skip_on_cran()
  expect_that(pcmsn <- bruvo.msn(pc, replen = rep(1, 10)), not(throws_error()))
  expect_equivalent(sort(unique(igraph::V(pcmsn$graph)$label)), sort(mll.levels(pc)))
  expect_that(plot_poppr_msn(pc, pcmsn), not(throws_error()))
  expect_that(plot_poppr_msn(pc, pcmsn, mlg = TRUE), not(throws_error()))
})