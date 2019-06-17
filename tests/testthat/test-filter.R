context("Multilocus Genotype Filter Tests")

lu <- function(x) length(unique(x))
grid_example <- matrix(c(1, 1, 5, 9, 9, 
                         4, 1, 1, 1, 4), 
                       ncol = 2)
rownames(grid_example) <- LETTERS[1:5]
colnames(grid_example) <- c("x", "y")
x  <- as.genclone(df2genind(grid_example, ploidy = 1))
xd <- dist(grid_example)
xdm <- as.matrix(xd)

set.seed(999)
gc <- as.snpclone(glSim(100, 0, n.snp.struc = 1e3, ploidy = 2, n.cores = 1L, parallel = FALSE))


test_that("multilocus genotype filtering algorithms work", {
  skip_on_cran()
  expect_equal(nmll(x), nInd(x))
  expect_equal(nmll(x), 5)
  
  # The following tests are from the paper
  fart <- mlg.filter(x, distance = xd, threshold = 4.51, algorithm = "f")
  aver <- mlg.filter(x, distance = xd, threshold = 4.51, algorithm = "a")
  near <- mlg.filter(x, distance = xd, threshold = 4.51, algorithm = "n")
  
  # Nearest should chain everything
  expect_equal(lu(near), 1)
  expect_equal(near, c(1L, 1L, 1L, 1L, 1L))
  
  # Average should make two groups
  expect_equal(lu(aver), 2)
  expect_equal(aver, c(4L, 4L, 1L, 1L, 1L))
  
  # Farthest should make three
  expect_equal(lu(fart), 3)
  expect_equal(fart, c(4L, 4L, 3L, 1L, 1L))
  
})

test_that("filter_stats works for genind", {
	skip_on_cran()
	res <- filter_stats(x, distance = xd, threshold = 100L, plot = TRUE, nclone = 2)
	maxres <- nmll(x) - 1
	
	expect_is(res, "list")
	expect_equal(names(res), c("farthest", "average", "nearest"))
	res <- lapply(res, "[[", "THRESHOLDS")
	expect_true(all(unlist(lapply(res, "[", 1:2), use.names = FALSE) == 3))
	expect_true(all(vapply(res, length, numeric(1)) == maxres))
	expect_lt(res$nearest[maxres], res$average[maxres])
	expect_lt(res$average[maxres], res$farthest[maxres])

  cpfart <- cutoff_predictor(res$farthest, 0.75)
  cpaver <- cutoff_predictor(res$average, 0.75)
  cpnear <- cutoff_predictor(res$nearest, 0.75)

  expect_equivalent(cpfart, 4)
  expect_equivalent(cpaver, 3.75)
  expect_equivalent(cpnear, 3.5)
})

test_that("filter_stats still plots if thresholds is selected", {
  skip_on_cran()
  tmp <- tempfile(pattern = "hay", fileext = ".pdf")
  expect_true(is.na(file.size(tmp)))
  pdf(tmp)
  filter_stats(x, distance = xd, threshold = 100L, plot = TRUE, nclone = 2, stats = "threshold")
  dev.off()
  expect_true(file.exists(tmp))
  expect_true(file.size(tmp) > 0)
  file.remove(tmp)
})

test_that("filter_stats works for snpclone", {
	skip_on_cran()
	res <- filter_stats(gc, distance = bitwise.dist, threshold = 100L, plot = TRUE)
	maxres <- nmll(gc) - 1
	
	expect_is(res, "list")
	expect_equal(names(res), c("farthest", "average", "nearest"))
	expect_is(res$nearest, "list")
	res <- lapply(res, "[[", "THRESHOLDS")
	expect_true(all(vapply(res, length, numeric(1)) == maxres))
	expect_lt(res$nearest[maxres], res$average[maxres])
	expect_lt(res$average[maxres], res$farthest[maxres])
})


test_that("mlg.filter can remember things", {
  skip_on_cran()
  
  # Default printing does not show code
  tda <- "(\\[t\\]).+?(\\[d\\]).+?(\\[a\\])"
  expect_true(all( !grepl(tda, capture.output(x)) ))
  expect_true(all( !grepl(tda, capture.output(gc)) ))
  suppressWarnings(mlg.filter(x) <- 0)
  suppressWarnings(mlg.filter(gc) <- 0)
  
  # Code is shown after modification and values are set to 0, nei and farthest
  expect_output(show(x), tda)
  original_vals <- "(0).+?(diss.dist).+?(farthest).+?"
  gc_original_vals <- "(0).+?(bitwise.dist).+?(farthest).+?"
  expect_output(show(x), original_vals)
  expect_output(print(x), original_vals)
  expect_output(show(gc), gc_original_vals)
  expect_output(print(gc), gc_original_vals)
})
  
  
test_that("supplied distance matrices work", {
  assign("x20150702210257_distance", xd)
  assign("x20150703173505_distance", bitwise.dist(gc, differences_only = TRUE))
  mlg.filter(x, distance = x20150702210257_distance) <- 0
  mlg.filter(gc, distance = x20150703173505_distance) <- 0
  expect_output(show(x), "x20150702210257_distance")
  expect_output(show(gc), "x20150703173505_distance")
  
  # choosing different algorithms work
  mlg.filter(x, algo = "average") <- 0
  mlg.filter(gc, algo = "average") <- 0
  expect_output(show(x), "average")
  expect_output(show(gc), "average")
  
  # and the distance is still persistent
  expect_output(show(x), "x20150702210257_distance")
  expect_output(show(gc), "x20150703173505_distance")

  # choosing different distances with parameters works
  mlg.filter(x, distance = diss.dist, percent = TRUE) <- 0
  expect_equal(nmll(x), 5)
  expect_output(show(x), "diss.dist")
  
  mlg.filter(gc, distance = bitwise.dist, percent = FALSE) <- 0
  expect_equal(nmll(gc), 100)
  expect_output(show(gc), "bitwise.dist")
  
  # algorithms are persistant
  expect_output(show(x), "average")
  expect_output(show(gc), "average")
  
  # Parameters are kept (default for diss.dist is percent = FALSE)
  mlg.filter(x) <- 0.75
  expect_equal(nmll(x), 3)
  expect_output(show(x), "0.75")
  
  # (default for bitwise.dist is percent = TRUE)
  mlg.filter(gc) <- 500
  expect_equal(nmll(gc), 81)
  expect_output(show(gc), "500")
  
  mlg.filter(x, distance = x20150702210257_distance) <- 4.51
  expect_equal(nmll(x), 2)
  
  mlg.filter(gc, distance = x20150703173505_distance) <- 0.25
  expect_equal(nmll(gc), 100)
  mlg.filter(gc) <- 0.45
  expect_lt(nmll(gc), 100)
  
  # An error is thrown if the distance is removed
  rm("x20150702210257_distance")
  rm("x20150703173505_distance")
  expect_error(mlg.filter(x) <- 0)
  expect_error(mlg.filter(gc) <- 0)
  
})


test_that("users can change parameter arguments", {

  skip_on_cran()
  mlg.filter(gc, distance = bitwise.dist, percent = TRUE) <- 0.25
  nmlg_gc_perc <- nmll(gc)
  expect_lt(nmlg_gc_perc, 100)
  expect_identical(distargs(gc@mlg), list(percent = TRUE))

  mlg.filter(gc, distance = bitwise.dist, percent = FALSE, euclidean = TRUE) <- 25
  nmlg_gc_eucl <- nmll(gc)
  # The number of mlg is less with euclidean than percentage
  expect_lt(nmlg_gc_eucl, nmlg_gc_perc)
  expect_identical(distargs(gc@mlg), list(percent = FALSE, euclidean = TRUE))

})

test_that("filtering algorithms imply diss.dist by default", {
  skip_on_cran()
  x <- structure(c(4, 1, 2, 2, 1, NA, 4, 4, 1, 1, 3, NA, 2, 1, 4, 2, 
  1, 4, 4, 2, 1, NA, 4, 4, 4, 1, 4, 2, 1, 4, 4, 2, 2, 1, 4, 4, 
  2, 1, 3, 2, 1, 4, 4, 4, 2, 1, 4, 4, 2, 1, 2, 2, 1, 4, 4, 3, 4, 
  1, 4, 4), .Dim = c(12L, 5L))
  y <- as.genclone(df2genind(x, ploidy = 1))
  mlg.filter(y) <- 0.001
  expect_gt(nmll(y, "original"), nmll(y, "contracted"))
  expect_output(show(y), "diss.dist")
})


context("mlg.filter error messages")

data("monpop", package = "poppr")
assign("x20160810_mon20", monpop[1:20])
assign("x20160810_let", matrix(letters[1:9], 3, 3))
assign("x20160810_neifun", function(x) nei.dist(genind2genpop(x, quiet = TRUE)))


test_that("internal filtering will throw an error for negative distances", {
  skip_on_cran()
  xdn <- xd
  xdn[1] <- -3
  expect_error(mlg.filter(x, distance = xdn, threshold = 4.51), "Distance matrix must not contain negative distances")
  xdn[1] <- -0.3
  expect_warning(.Call("neighbor_clustering", as.matrix(xdn), mll(x), 4.51, "f", 1L), "The data resulted in a negative or invalid distance or cluster id")
})

test_that("a warning is thrown if the user specifies more than one thread.", {
  skip_on_cran()
  expect_warning(mlg.filter(x, distance = xd, threshold = 4.51, threads = 2L))
})

test_that("Infinite distances will produce an error", {
  skip_on_cran()
  xdn <- xd
  xdn[1] <- Inf
  expect_error(mlg.filter(x, distance = xdn, threshold = 4.51), "Data set contains missing or invalid distances")
})


test_that("mlg.filter errors when distance matrix is not square", {
  skip_on_cran()
  mat <- matrix(runif(nInd(monpop)*10), nInd(monpop), 10)
  msg <- "must be a square matrix"
  expect_error(mlg.filter(monpop, distance = mat) <- 0, msg)
})

test_that("mlg.filter errors when distance matrix contains missing data", {
  skip_on_cran()
  expect_error(mlg.filter(x20160810_mon20, distance = nei.dist) <- 0, "contains missing data")
})

test_that("mlg.filter errors when the distance matrix is not of equal size", {
  skip_on_cran()
  expect_error(mlg.filter(monpop, distance = diss.dist(x20160810_mon20)) <- 0, "distance matrix \\(20\\)")
  expect_error(mlg.filter(monpop, distance = x20160810_neifun) <- 0, "samples, not populations")
})

test_that("mlg.filter throws a warning when attempting to use 'mean' with diss.dist", {
  skip_on_cran()
  expect_warning(mlg.filter(x20160810_mon20, distance = diss.dist, missing = "mean") <- 0)
})

test_that("mlg.filter throws an error if the threshold is not a number", {
  skip_on_cran()
  expect_error(mlg.filter(x) <- "A", "Threshold must be")
})

test_that("mlg.filter throws an error if the distance is not numeric", {
  skip_on_cran()
  expect_error(mlg.filter(x, distance = xdm > 4) <- 1, "Distance matrix must be")
})

rm("x20160810_mon20")
rm("x20160810_let")
rm("x20160810_neifun")
