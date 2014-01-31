context("Analytical value tests")

test_that("Bruvo's distance works as expected.", {
  testdf  <- data.frame(test = c("00/20/23/24", "20/24/26/43"))
  testgid <- df2genind(testdf, ploidy = 4, sep = "/")
  addloss <- as.vector(bruvo.dist(testgid, add = FALSE, loss = FALSE))
  ADDloss <- as.vector(bruvo.dist(testgid, add = TRUE, loss = FALSE))
  addLOSS <- as.vector(bruvo.dist(testgid, add = FALSE, loss = TRUE))
  ADDLOSS <- as.vector(bruvo.dist(testgid, add = TRUE, loss = TRUE))
  expect_that(addloss, equals(0.46875000000000))
  expect_that(ADDloss, equals(0.34374987334013))
  expect_that(addLOSS, equals(0.458333164453506))
  expect_that(ADDLOSS, equals(0.401041518896818))
})

test_that("Dissimilarity distance works as expected.", {
  data(nancycats)
  nan1 <- popsub(nancycats, 1)
  nanmat <- diss.dist(nan1, mat = TRUE)
  expect_that(diss.dist(nan1), is_a("dist"))
  expect_that(nanmat, is_a("matrix"))
  expect_that(diss.dist(nan1, mat = TRUE, frac = FALSE), equals(nanmat * 2 * 9))
  expect_that(nanmat[2, 1], equals(0.222222222222222))
})

test_that("Index of association works as expected.", {
  data(Aeut)
  res <- c(Ia = 14.3707995986407, rbarD = 0.270617053778004)
  expect_that(ia(Aeut), equals(res))
})