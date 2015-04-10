context("Analytical value tests")

test_that("Bruvo's distance works as expected.", {
  testdf  <- data.frame(test = c("00/20/23/24", "20/24/26/43"))
  testgid <- df2genind(testdf, ploidy = 4, sep = "/")
  addloss <- as.vector(bruvo.dist(testgid, add = FALSE, loss = FALSE))
  ADDloss <- as.vector(bruvo.dist(testgid, add = TRUE, loss = FALSE))
  addLOSS <- as.vector(bruvo.dist(testgid, add = FALSE, loss = TRUE))
  ADDLOSS <- as.vector(bruvo.dist(testgid, add = TRUE, loss = TRUE))
  expect_that(addloss, equals(0.46875000000000))
  expect_that(addLOSS, equals(0.34374987334013))
  expect_that(ADDloss, equals(0.458333164453506))
  expect_that(ADDLOSS, equals(0.401041518896818))
})

test_that("Infinite Alleles Model works.",{
  x <- structure(list(V3 = c("228/236/242", "000/211/226"), 
                      V6 = c("190/210/214", "000/190/203")), 
                 .Names = c("V3", "V6"), row.names = c("1", "2"), 
                 class = "data.frame")
  gid <- df2genind(x, sep = "/", ploidy = 3)
  res <- bruvo.dist(gid, replen = c(2/10000,2/10000), add = FALSE, loss = FALSE)
  expect_that(as.numeric(res), equals(0.833333333333333))
})

test_that("Dissimilarity distance works as expected.", {
  data(nancycats, package = "adegenet")
  nan1 <- popsub(nancycats, 1)
  nanmat <- diss.dist(nan1, mat = TRUE)
  expect_that(diss.dist(nan1), is_a("dist"))
  expect_that(nanmat, is_a("matrix"))
  #expect_that(diss.dist(nan1, mat = TRUE, percent = FALSE), equals(nanmat*2*9))
  #expect_that(nanmat[2, 1], equals(0.222222222222222))
  expect_that(diss.dist(nan1, mat = TRUE, percent = TRUE), equals((nanmat/2)/9))
  expect_that(nanmat[2, 1], equals(4))
})

test_that("Index of association works as expected.", {
  data(Aeut, package = "poppr")
  res <- c(Ia = 14.3707995986407, rbarD = 0.270617053778004)
  expect_that(ia(Aeut), equals(res))
})

test_that("Internal function fix_negative_branch works as expected.", {
  the_tree <- structure(list(edge = structure(c(9L, 10L, 11L, 12L, 13L, 14L,
      14L, 13L, 12L, 11L, 10L, 9L, 9L, 10L, 11L, 12L, 13L, 14L, 2L, 3L, 6L, 7L,
      8L, 1L, 5L, 4L), .Dim = c(13L, 2L)), 
      edge.length = c(0, 0.0625, 0.0625, 0.09375, 0.15, -0.25, 0.25, -0.15,
      -0.09375, -0.0625, 0, 0, 0),
      tip.label = c("2340_50156.ab1 ", "2340_50149.ab1 ", "2340_50674.ab1 ",
      "2370_45312.ab1 ", "2340_50406.ab1 ", "2370_45424.ab1 ", "2370_45311.ab1
      ", "2370_45521.ab1 "), 
      Nnode = 6L
    ),
    .Names = c("edge", "edge.length", "tip.label", "Nnode"), 
    class = "phylo",
    order = "cladewise")
  fix_tree <- poppr:::fix_negative_branch(the_tree)
  # Not all branch lengths are positive
  expect_false(min(the_tree$edge.length) >= 0)
  # After fix, all branch lengths are positive
  expect_true(min(fix_tree$edge.length) >= 0)
  # The difference from fixed and unfixed is unfixed. This indicates that the
  # clones were set to zero and the fix set the branch lengths in the correct 
  # order.
  expect_equivalent(min(fix_tree$edge.length - the_tree$edge.length), min(the_tree$edge.length))
})