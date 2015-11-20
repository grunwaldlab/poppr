context("Minimum spanning network functions")

options(warn = -1)
ucl <- function(x){
  unclass(x$graph)[-10]
}

set.seed(9005)

gend <-new("genind"
           , tab = structure(c(0.5, 0.5, 0, 0, 0.5, 0.5, 1, 1, 0, 0.5, 0.5, 0, 1,
                               0.5, 0.5, 1)*2, .Dim = c(4L, 4L), .Dimnames = list(c("1", "2",
                                                                                  "3", "4"), c("loc-1.1", "loc-1.2", "loc-2.1", "loc-2.2")))
           , loc.names = structure(c("loc-1", "loc-2"), .Names = c("loc-1", "loc-2"))
           , loc.fac = structure(c(1L, 1L, 2L, 2L), .Label = c("loc-1", "loc-2"), class = "factor")
           , loc.nall = structure(c(2L, 2L), .Names = c("loc-1", "loc-2"))
           , all.names = structure(list(`loc-1` = structure(c("3", "4"), .Names = c("1", "2"
           )), `loc-2` = structure(c("3", "4"), .Names = c("1", "2"))), .Names = c("loc-1",
                                                                              "loc-2"))
           , call = NULL
           , ind.names = structure(c("", "", "", ""), .Names = c("1", "2", "3", "4"))
           , pop = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor")
           , pop.names = structure(c("1", "2", "3", "4"))
           , ploidy = 2L
           , type = "codom"
           , other = NULL
)

set.seed(9005)

gend_single <- new("genind"
    , tab = structure(c(0.5, 0.5, 0, 0, 0.5, 0.5, 1, 1, 0, 0.5, 0.5, 0, 1, 0.5, 0.5, 1)*2, .Dim = c(4L, 4L), .Dimnames = list(c("1", "2",
"3", "4"), c("loc-1.1", "loc-1.2", "loc-2.1", "loc-2.2")))
    , loc.names = structure(c("loc-1", "loc-2"), .Names = c("loc-1", "loc-2"))
    , loc.fac = structure(c(1L, 1L, 2L, 2L), .Label = c("loc-1", "loc-2"), class = "factor")
    , loc.nall = structure(c(2L, 2L), .Names = c("loc-1", "loc-2"))
    , all.names = structure(list(`loc-1` = structure(c("3", "4"), .Names = c("1", "2")), `loc-2` = structure(c("3", "4"), .Names = c("1", "2"))), .Names = c("loc-1",
"loc-2"))
    , call = NULL
    , ind.names = structure(c("", "", "", ""), .Names = c("1", "2", "3", "4"))
    , pop = structure(c(1,1,1,1), .Label = c("1"), class = "factor")
    , pop.names = structure(c("1","1","1","1"))
    , ploidy = 2L
    , type = "codom"
    , other = NULL
)

no_ties <- structure(list(graph = structure(list(4, FALSE, c(3, 2, 3), c(0,
1, 2), c(1, 0, 2), c(0, 1, 2), c(0, 0, 0, 1, 3), c(0, 1, 2, 3,
3), list(c(1, 0, 1), structure(list(), .Names = character(0)),
    structure(list(name = c("1", "2", "3", "4"), size = c(1L,
    1L, 1L, 1L), shape = c("pie", "pie", "pie", "pie"), pie = list(
        structure(1L, .Names = "1"), structure(1L, .Names = "2"),
        structure(1L, .Names = "3"), structure(1L, .Names = "4")),
        pie.color = list(structure("#4C00FFFF", .Names = "1"), structure("#00E5FFFF", .Names = "2"), structure("#00FF4DFF", .Names = "3"),
            structure("#FFFF00FF", .Names = "4")), label = c("MLG.3", "MLG.4", "MLG.2",
        "MLG.1")), .Names = c("name", "size", "shape", "pie",
    "pie.color", "label")), structure(list(weight = c(0.125,
    0.125, 0.125), color = c("#434343", "#434343", "#434343"),
        width = c(8, 8, 8)), .Names = c("weight", "color", "width"
    )))), class = "igraph"), populations = structure(c("1", "2",
"3", "4")), colors = structure(c("#4C00FFFF",
"#00E5FFFF", "#00FF4DFF", "#FFFF00FF"), .Names = c("1","2","3","4"))), .Names = c("graph",
"populations", "colors"))

no_ties_single <- structure(list(graph = structure(list(4, FALSE, c(3, 2, 3), c(0,
1, 2), c(1, 0, 2), c(0, 1, 2), c(0, 0, 0, 1, 3), c(0, 1, 2, 3,
3), list(c(1, 0, 1), structure(list(), .Names = character(0)),
    structure(list(name = c("1", "2", "3", "4"), size = c(1L,
    1L, 1L, 1L), color = c("#4C00FFFF", "#4C00FFFF", "#4C00FFFF",
    "#4C00FFFF"), label = c("MLG.3", "MLG.4", "MLG.2", "MLG.1"
    )), .Names = c("name", "size", "color", "label")), structure(list(
        weight = c(0.125, 0.125, 0.125), color = c("#434343",
        "#434343", "#434343"), width = c(8, 8, 8)), .Names = c("weight",
    "color", "width")))), class = "igraph"), populations = "1",
    colors = "#4C00FFFF"), .Names = c("graph", "populations",
"colors"))

ties <- structure(list(graph = structure(list(4, FALSE, c(3, 2, 3, 1),
    c(0, 1, 2, 0), c(3, 1, 0, 2), c(3, 0, 1, 2), c(0, 0, 1, 2,
    4), c(0, 2, 3, 4, 4), list(c(1, 0, 1), structure(list(), .Names = character(0)),
        structure(list(name = c("1", "2", "3", "4"), size = c(1L,
        1L, 1L, 1L), shape = c("pie", "pie", "pie", "pie"), pie = list(
            structure(1L, .Names = "1"), structure(1L, .Names = "2"),
            structure(1L, .Names = "3"), structure(1L, .Names = "4")),
            pie.color = list(structure("#4C00FFFF", .Names = "1"), structure("#00E5FFFF", .Names = "2"),
                structure("#00FF4DFF", .Names = "3"), structure("#FFFF00FF", .Names = "4")),
            label = c("MLG.3", "MLG.4", "MLG.2",
            "MLG.1")), .Names = c("name", "size", "shape", "pie",
        "pie.color", "label")), structure(list(weight = c(0.125,
        0.125, 0.125, 0.125), color = c("#434343", "#434343",
        "#434343", "#434343"), width = c(8, 8, 8, 8)), .Names = c("weight",
        "color", "width")))), class = "igraph"), populations = structure(c("1",
"2", "3", "4")), colors = structure(c("#4C00FFFF",
"#00E5FFFF", "#00FF4DFF", "#FFFF00FF"), .Names = c("1","2","3","4"))), .Names = c("graph",
"populations", "colors"))

ties_single <- structure(list(graph = structure(list(4, FALSE, c(3, 2, 3, 1),
    c(0, 1, 2, 0), c(3, 1, 0, 2), c(3, 0, 1, 2), c(0, 0, 1, 2,
    4), c(0, 2, 3, 4, 4), list(c(1, 0, 1), structure(list(), .Names = character(0)),
        structure(list(name = c("1", "2", "3", "4"), size = c(1L,
        1L, 1L, 1L), color = c("#4C00FFFF", "#4C00FFFF", "#4C00FFFF",
        "#4C00FFFF"), label = c("MLG.3", "MLG.4", "MLG.2", "MLG.1"
        )), .Names = c("name", "size", "color", "label")), structure(list(
            weight = c(0.125, 0.125, 0.125, 0.125), color = c("#434343",
            "#434343", "#434343", "#434343"), width = c(8, 8,
            8, 8)), .Names = c("weight", "color", "width")))), class = "igraph"),
    populations = "1", colors = "#4C00FFFF"), .Names = c("graph",
"populations", "colors"))

if (packageVersion("igraph") >= package_version("1.0.0")){
    no_ties$graph <- igraph::upgrade_graph(no_ties$graph)
    no_ties_single$graph <- igraph::upgrade_graph(no_ties_single$graph)
    ties$graph <- igraph::upgrade_graph(ties$graph)
    ties_single$graph <- igraph::upgrade_graph(ties_single$graph)
}

options(warn = 0)

test_that("Minimum spanning networks can properly account for tied edges", {

  # Test Bruvo.msn
  set.seed(9005)
  expect_equal(ucl(bruvo.msn(gend, replen=c(1,1))), ucl(no_ties))
  set.seed(9005)
  expect_equal(ucl(bruvo.msn(gend, replen=c(1,1), include.ties = TRUE)), ucl(ties))

  # Test poppr.msn
  set.seed(9005)
  expect_equal(ucl(poppr.msn(gend, distmat=bruvo.dist(gend,replen=c(1,1)))), ucl(no_ties))
  set.seed(9005)
  expect_equal(ucl(poppr.msn(gend, distmat=bruvo.dist(gend,replen=c(1,1)), include.ties = TRUE)), ucl(ties))

  # Test both for single populations sets
  set.seed(9005)
  expect_equal(ucl(bruvo.msn(gend_single, replen=c(1,1))), ucl(no_ties_single))
  set.seed(9005)
  expect_equal(ucl(bruvo.msn(gend_single, replen=c(1,1), include.ties = TRUE)), ucl(ties_single))
  set.seed(9005)
  expect_equal(ucl(poppr.msn(gend_single, distmat=bruvo.dist(gend_single,replen=c(1,1)))), ucl(no_ties_single))
  set.seed(9005)
  expect_equal(ucl(poppr.msn(gend_single, distmat=bruvo.dist(gend_single,replen=c(1,1)), include.ties = TRUE)), ucl(ties_single))

})

test_that("Minimum spanning networks also collapse MLGs", {
  skip_on_cran()
  gend        <- as.genclone(gend)
  gend_single <- as.genclone(gend_single)
  gend_bruvo  <- bruvo.dist(gend, replen = c(1, 1))

  gmsnt  <- bruvo.msn(gend, replen = c(1, 1), threshold = 0.15)
  pgmsnt <- poppr.msn(gend, distmat = gend_bruvo, threshold = 0.15)

  expect_identical(igraph::V(gmsnt$graph)$pie, igraph::V(pgmsnt$graph)$pie)
  expect_identical(igraph::V(gmsnt$graph)$name, igraph::V(pgmsnt$graph)$name)
  expect_identical(igraph::E(gmsnt$graph)$weight, igraph::E(pgmsnt$graph)$weight)
  
  
  gmsn   <- bruvo.msn(gend, replen = c(1, 1), showplot = FALSE)
  pgmsn  <- poppr.msn(gend, distmat = gend_bruvo, showplot = FALSE)

  expect_identical(igraph::V(gmsn$graph)$pie, igraph::V(pgmsn$graph)$pie)
  expect_identical(igraph::V(gmsn$graph)$name, igraph::V(pgmsn$graph)$name)
  expect_identical(igraph::E(gmsn$graph)$weight, igraph::E(pgmsn$graph)$weight)

  expect_equal(length(igraph::V(gmsnt$graph)), 2)
  expect_equal(length(igraph::V(gmsn$graph)), 4)

  sgmsnt  <- bruvo.msn(gend_single, replen = c(1, 1), threshold = 0.15)
  psgmsnt <- poppr.msn(gend_single, distmat = gend_bruvo, threshold = 0.15)

  expect_identical(igraph::V(sgmsnt$graph)$pie, igraph::V(psgmsnt$graph)$pie)
  expect_identical(igraph::V(sgmsnt$graph)$name, igraph::V(psgmsnt$graph)$name)
  expect_identical(igraph::E(sgmsnt$graph)$weight, igraph::E(psgmsnt$graph)$weight)
  
  sgmsn   <- bruvo.msn(gend_single, replen = c(1, 1), showplot = FALSE)
  psgmsn  <- poppr.msn(gend_single, distmat = gend_bruvo, showplot = FALSE)

  expect_identical(igraph::V(sgmsn$graph)$pie, igraph::V(psgmsn$graph)$pie)
  expect_identical(igraph::V(sgmsn$graph)$name, igraph::V(psgmsn$graph)$name)
  expect_identical(igraph::E(sgmsn$graph)$weight, igraph::E(psgmsn$graph)$weight)

  expect_equal(length(igraph::V(sgmsnt$graph)), 2)
  expect_equal(length(igraph::V(sgmsn$graph)), 4)

  expect_output(plot_poppr_msn(gend, gmsnt, palette = "cm.colors"), "")
  expect_output(plot_poppr_msn(gend_single, sgmsnt, palette = "cm.colors"), "")
})


test_that("Filtered minimum spanning networks retain original names", {
  skip_on_cran()
  # setup ----------------------------------------------------
  grid_example <- matrix(c(1, 4,
                           1, 1,
                           5, 1,
                           9, 1,
                           9, 4), 
                         ncol = 2,
                         byrow = TRUE)
  rownames(grid_example) <- LETTERS[1:5]
  colnames(grid_example) <- c("x", "y")
  grid_new <- rbind(grid_example, 
                    new = c(5, NA), 
                    mut = c(5, 2)
                    )
  x <- as.genclone(df2genind(grid_new, ploidy = 1))
  indNames(x)
  ## [1] "A"   "B"   "C"   "D"   "E"   "new" "mut"
  
  raw_dist <- function(x){
    dist(genind2df(x, usepop = FALSE))
  }
  (xdis <- raw_dist(x))
  
  # normal ---------------------------------------------------
  set.seed(9001)
  g1 <- poppr.msn(x, xdis, include.ties = TRUE, showplot = FALSE,
                  vertex.label.color = "firebrick", vertex.label.font = 2)
  all_names <- igraph::V(g1$graph)$name
  ## [1] "A"   "B"   "C"   "D"   "E"   "new" "mut"
  
  # filtered ---------------------------------------------------
  set.seed(9001)
  g1.1 <- poppr.msn(x, xdis, threshold = 1, include.ties = TRUE, showplot = FALSE,
                  vertex.label.color = "firebrick", vertex.label.font = 2)
  cc_names <- igraph::V(g1.1$graph)$name
  ## [1] "A"   "B"   "C"   "D"   "E"   "new"
  expect_identical(cc_names, head(all_names, -1))
})

data("partial_clone")
pc <- as.genclone(partial_clone)

test_that("Minimum spanning networks can subset populations", {
  bpc  <- bruvo.dist(pc, replen = rep(1, 10))
  bmsn <- bruvo.msn(pc, replen = rep(1, 10), showplot = FALSE)
  pmsn <- poppr.msn(pc, bpc, showplot = FALSE)
  expect_identical(ucl(bmsn), ucl(pmsn))
  
  bmsn12 <- bruvo.msn(pc, replen = rep(1, 10), sublist = 1:2, showplot = FALSE)
  pmsn12 <- poppr.msn(pc, bpc, sublist = 1:2, showplot = FALSE)
  expect_identical(ucl(bmsn12), ucl(pmsn12))

  bmsn1 <- bruvo.msn(pc, replen = rep(1, 10), sublist = 1, showplot = FALSE)
  pmsn1 <- poppr.msn(pc, bpc, sublist = 1, showplot = FALSE)
  expect_identical(ucl(bmsn1), ucl(pmsn1))
  
})


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

test_that("vectors can be used to color graphs", {
  skip_on_cran()
  data(Aeut)
  A.dist <- diss.dist(Aeut)
  
  # Graph it.
  A.msn <- poppr.msn(Aeut, A.dist, gadj=15, vertex.label=NA, showplot = FALSE)
  unpal <- c("black", "orange")
  fpal  <- function(x) unpal
  npal  <- setNames(unpal, c("Athena", "Mt. Vernon"))
  xpal  <- c(npal, JoMo = "awesome")
  # Using palette without names
  uname_pal  <- plot_poppr_msn(Aeut, A.msn, palette = unpal)$colors
  # Using palette with function
  fun_pal    <- plot_poppr_msn(Aeut, A.msn, palette = fpal)$colors
  # Using palette with names
  name_pal   <- plot_poppr_msn(Aeut, A.msn, palette = npal[2:1])$colors
  # Using palette with extra names
  xname_pal  <- plot_poppr_msn(Aeut, A.msn, palette = xpal)$colors
  
  expect_identical(uname_pal, npal)
  expect_identical(fun_pal, npal)
  expect_identical(name_pal, npal)
  expect_identical(xname_pal, npal)
})
