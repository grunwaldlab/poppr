context("win.ia tests")
set.seed(999)
chrom_pos <- vapply(1:10, function(i) sort(sample(1e3, 100)), integer(100)) %>%
  as.vector()
chromo    <- rep(1:10, each = 100)
set.seed(999)
x <- glSim(n.ind = 10, n.snp.nonstruc = 5e2, n.snp.struc = 5e2, ploidy = 2, 
           parallel = FALSE)


test_that("win.ia will throw an error if duplicate positions are found", {
  options(poppr.debug = TRUE)
  expect_output(x.naive     <- win.ia(x, name_window = TRUE), "[|=]{2,}")
  expect_equal(length(x.naive), 10L)
  expect_named(x.naive, as.character(100 * (1:10)))
  expect_null(names(win.ia(x, quiet = TRUE, name_window = FALSE)))
  position(x) <- chrom_pos
  expect_error(win.ia(x), "chromosome")
  options(poppr.debug = FALSE)
})  

test_that("win.ia will throw a warning if chromosome_buffer is specified", {
  expect_warning(x.naive <- win.ia(x, quiet = TRUE, chromosome_buffer = FALSE), "deprecated")
})

test_that("win.ia will use chromosome structure", {
  skip_on_cran()
  position(x)   <- chrom_pos
  chromosome(x) <- chromo
  x.chrom.bt    <- win.ia(x, quiet = TRUE)
  expect_equal(length(x.chrom.bt), 100L)
  winnames <- paste(rep(1:10, each = 10), rep(100*(1:10), 10), sep = ".")
  expect_equal(names(x.chrom.bt), winnames)
})
  
test_that("win.ia will always start at the beginning of the chromosome", {
  skip_on_cran()
  position(x)   <- chrom_pos
  chromosome(x) <- chromo
  x.chrom.bf    <- win.ia(x, window = 300L, quiet = TRUE)
  x.by.chrom    <- lapply(levels(chromosome(x)), function(i) win.ia(x[, chromosome(x) == i], window = 300L, quiet = TRUE))
  # There should be 4 windows per chromosome.
  expect_equal(length(x.chrom.bf), 40L)
  # Using the function with and without chromosome structure should not matter.
  expect_equal(x.chrom.bf, unlist(x.by.chrom, use.names = FALSE), check.attributes = FALSE)
  
})
