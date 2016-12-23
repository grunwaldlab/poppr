context("win.ia tests")
options(poppr.debug = FALSE)
set.seed(999)
chrom_pos <- vapply(1:10, function(i) sort(sample(1e3, 100)), integer(100)) %>%
  as.vector()
chromo    <- rep(1:10, each = 100)
real_pos  <- chrom_pos + 1000*(chromo - 1) + 1#+ 100*(chromo - 1) + 1  
set.seed(999)
x <- glSim(n.ind = 10, n.snp.nonstruc = 5e2, n.snp.struc = 5e2, ploidy = 2, 
           parallel = FALSE)


test_that("win.ia will throw an error if duplicate positions are found", {
  x.naive     <- win.ia(x)
  position(x) <- chrom_pos
  expect_equal(length(x.naive), 10L)
  expect_error(win.ia(x), "chromosome")
})  

test_that("win.ia will use chromosome structure", {
  skip_on_cran()
  position(x)   <- chrom_pos
  chromosome(x) <- chromo
  x.chrom.bt    <- win.ia(x)
  expect_equal(sum(is.na(x.chrom.bt) & !is.nan(x.chrom.bt)), 9L)
  expect_equal(length(x.chrom.bt), 109L)
})
  
test_that("win.ia will take into account chromosome buffer", {
  skip_on_cran()
  position(x)   <- chrom_pos
  chromosome(x) <- chromo
  x.chrom.bf    <- win.ia(x, chromosome_buffer = FALSE)
  chromosome(x) <- NULL
  position(x)   <- real_pos
  x.real        <- win.ia(x)
  expect_equal(length(x.chrom.bf), 100L)
  expect_equal(length(x.chrom.bf), length(x.real))
  expect_equal(x.chrom.bf, x.real, check.attributes = FALSE)
})

