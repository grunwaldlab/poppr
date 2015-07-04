context("Multilocus Genotype Filter Tests")

lu <- function(x) length(unique(x))
grid_example <- matrix(c(1, 1, 5, 9, 9, 
                         4, 1, 1, 1, 4), 
                       ncol = 2)
rownames(grid_example) <- LETTERS[1:5]
colnames(grid_example) <- c("x", "y")
x  <- as.genclone(df2genind(grid_example, ploidy = 1))
xd <- dist(grid_example)

set.seed(999)
gc <- as.snpclone(glSim(100, 0, n.snp.struc = 1e3, ploidy = 2, parallel = FALSE))


test_that("multilocus genotype filtering algorithms work", {
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

test_that("mlg.filter can remember things", {
  skip_on_cran()
  
  # Default printing does not show code
  tda <- "(\\[t\\]).+?(\\[d\\]).+?(\\[a\\])"
  expect_true(all( !grepl(tda, capture.output(x)) ))
  expect_true(all( !grepl(tda, capture.output(gc)) ))
  suppressWarnings(mlg.filter(x) <- 0)
  suppressWarnings(mlg.filter(gc) <- 0)
  
  # Code is shown after modification and values are set to 0, nei and farthest
  expect_output(x, tda)
  original_vals <- "(0).+?(nei.dist).+?(farthest).+?"
  gc_original_vals <- "(0).+?(bitwise.dist).+?(farthest).+?"
  expect_output(x, original_vals)
  expect_output(gc, gc_original_vals)
  
  # supplied distance matrices work
  assign("x20150702210257_distance", xd, envir = .GlobalEnv)
  assign("x20150703173505_distance", bitwise.dist(gc, differences_only = TRUE), envir = .GlobalEnv)
  mlg.filter(x, distance = x20150702210257_distance) <- 0
  mlg.filter(gc, distance = x20150703173505_distance) <- 0
  expect_output(x, "x20150702210257_distance")
  expect_output(gc, "x20150703173505_distance")
  
  # choosing different algorithms work
  mlg.filter(x, algo = "average") <- 0
  mlg.filter(gc, algo = "average") <- 0
  expect_output(x, "average")
  expect_output(gc, "average")
  
  # and the distance is still persistent
  expect_output(x, "x20150702210257_distance")
  expect_output(gc, "x20150703173505_distance")

  # choosing different distances with parameters works
  mlg.filter(x, distance = diss.dist, percent = TRUE) <- 0
  expect_equal(nmll(x), 5)
  expect_output(x, "diss.dist")
  
  mlg.filter(gc, distance = bitwise.dist, percent = FALSE) <- 0
  expect_equal(nmll(gc), 100)
  expect_output(gc, "bitwise.dist")
  
  # algorithms are persistant
  expect_output(x, "average")
  expect_output(gc, "average")
  
  # Parameters are kept (default for diss.dist is percent = FALSE)
  mlg.filter(x) <- 0.75
  expect_equal(nmll(x), 3)
  expect_output(x, "0.75")
  
  # (default for bitwise.dist is percent = TRUE)
  mlg.filter(gc) <- 500
  expect_equal(nmll(gc), 81)
  expect_output(gc, "500")
  
  mlg.filter(x, distance = x20150702210257_distance) <- 4.51
  expect_equal(nmll(x), 2)
  
  mlg.filter(gc, distance = x20150703173505_distance) <- 0.25
  expect_equal(nmll(gc), 100)
  mlg.filter(gc) <- 0.45
  expect_equal(nmll(gc), 58)
  
  # An error is thrown if the distance is removed
  rm("x20150702210257_distance", envir = .GlobalEnv)
  rm("x20150703173505_distance", envir = .GlobalEnv)
  expect_error(mlg.filter(x) <- 0)
  expect_error(mlg.filter(gc) <- 0)
  
})

