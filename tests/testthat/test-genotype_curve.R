context("Genotype Curve")
genos <- c("A/A", "A/B", "A/C", "B/B", "B/C", "C/C")

set.seed(500)
v <- sample(genos, 100, replace = TRUE)
w <- c(rep(genos[2], 99), genos[3])           # found by cutoff
x <- c(rep(genos[1], 98), genos[3], genos[2]) # found by MAF
y <- c(rep(genos[1], 99), genos[2])           # found by both
z <- sample(genos, 100, replace = TRUE)
suppressWarnings({
  dat <- df2genind(data.frame(v = v, w = w, x = x, y = y, z = z), sep = "/") %>% 
    as.genclone %>%
    clonecorrect(strata = NA)
})
dat2 <- dat
dat2@tab[25, 5:6] <- NA

test_that("genotype_curve's drop functionality works by default", {
  skip_on_cran()
  expect_message(x <- genotype_curve(dat[-nInd(dat)], plot = FALSE, quiet = TRUE), "Dropping monomorphic loci: w, y")
  expect_message(y <- genotype_curve(dat[-nInd(dat2)], plot = FALSE, quiet = TRUE), "Dropping monomorphic loci: w, y")
  expect_equal(ncol(x), 2L)
  expect_equal(ncol(y), 2L)
})

test_that("dropna can be switched off", {
  skip_on_cran()
  expect_message(x <- genotype_curve(dat2[-nInd(dat2)], dropna = FALSE, plot = FALSE, quiet = TRUE), "Dropping monomorphic loci: y")
  expect_equal(ncol(x), 3L)
})

test_that("genotype_curve's drop can be turned off", {
  skip_on_cran()
  expect_output(x <- genotype_curve(dat[-nInd(dat)],  drop = FALSE, plot = FALSE), "/4 loci")
  expect_output(y <- genotype_curve(dat[-nInd(dat2)], drop = FALSE, plot = FALSE), "/4 loci")
  expect_equal(ncol(x), 4L)
  expect_equal(ncol(y), 4L)
})
