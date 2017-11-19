
ucl <- function(x){
  unclass(x$graph)[-10]
}

set.seed(9005)

gend <- new("genind"
           , tab = structure(c(1L, 1L, 0L, 0L, 
                               1L, 1L, 2L, 2L,
                               0L, 1L, 1L, 0L, 
                               2L, 1L, 1L, 2L), 
              .Dim = c(4L, 4L), 
              .Dimnames = list(c("1", "2", "3", "4"), 
                               c("loc-1.1", "loc-1.2", "loc-2.1", "loc-2.2")))
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

gend_gc          <- as.genclone(gend)
gend_bruvo       <- bruvo.dist(gend, replen = c(1, 1))
gend_single      <- gend
pop(gend_single) <- rep(1, 4)

#' Identical population and color fields
#'
#' @param g1 graph 1
#' @param g2 graph 2
#'
#' @noRd
expect_identical_metadata <- function(g1, g2){
  eval(bquote(expect_identical(.(g1[-1]), .(g2[-1]))))
}

#' Vertex counts
#'
#' @param g a graph
#' @param n an integer equal to the number of vertices in the graph
#'
#' @noRd
expect_vcount <- function(g, n){
  eval(bquote(expect_equal(igraph::vcount(.(g$graph)), .(n))))
}

#' Edge counts
#'
#' @param g a graph
#' @param n an integer equal to the number of edges in the graph
#'
#' @noRd
expect_ecount <- function(g, n){
  eval(bquote(expect_equal(igraph::ecount(.(g$graph)), .(n))))
}

#' Identical attributes
#'
#' @param g1 graph 1
#' @param g2 graph 2
#' @param a  a character specifying which attribute to test
#'
#' @noRd
expect_identical_vertex_attr <- function(g1, g2, a){
  eval(bquote(expect_identical(igraph::vertex_attr(.(g1)$graph, .(a)), 
                               igraph::vertex_attr(.(g2)$graph, .(a))
  )
  ))
}

#' Identical attributes
#'
#' @param g1 graph 1
#' @param g2 graph 2
#' @param a  a character specifying which attribute to test
#'
#' @noRd
expect_identical_edge_attr <- function(g1, g2, a){
  eval(bquote(expect_identical(igraph::edge_attr(.(g1)$graph, .(a)), 
                               igraph::edge_attr(.(g2)$graph, .(a))
  )
  ))
}

#' Attribute names
#'
#' @param g a graph
#' @param a a character vector specifying all the named attributes
#'
#' @noRd
expect_vertex_attr <- function(g, a){
  eval(bquote(expect_equal(igraph::vertex_attr_names(.(g$graph)), .(a))))
}

#' Testing distance tables
#'
#' @param g a graph
#' @param d a tabulaton of the distances between nodes. For example, a graph
#'   with four nodes connected in a line would have 3:1, indicating that there
#'   are three pairs that can be reached with one step, two with two steps, and
#'   one pair needs three steps.
#'
#' @noRd
expect_distance_table <- function(g, d){
  eval(bquote(expect_equal(igraph::distance_table(.(g$graph))$res, .(d))))
}


context("Tied MSN edge tests")

bruv_no_ties  <- bruvo.msn(gend, replen = c(1,1), showplot = FALSE)
bruv_ties     <- bruvo.msn(gend, replen = c(1,1), showplot = FALSE, include.ties = TRUE)
pmsn_no_ties  <- poppr.msn(gend, distmat = gend_bruvo, showplot = FALSE)
pmsn_ties     <- poppr.msn(gend, distmat = gend_bruvo, showplot = FALSE, include.ties = TRUE)
sbruv_no_ties <- bruvo.msn(gend_single, replen = c(1,1), showplot = FALSE)
sbruv_ties    <- bruvo.msn(gend_single, replen = c(1,1), showplot = FALSE, include.ties = TRUE)
spmsn_no_ties <- poppr.msn(gend_single, distmat = gend_bruvo, showplot = FALSE)
spmsn_ties    <- poppr.msn(gend_single, distmat = gend_bruvo, showplot = FALSE, include.ties = TRUE)
pienames      <- c("name", "size", "shape", "pie", "pie.color", "label")
nopienames    <- c("name", "size", "color", "label")

test_that("bruvo.msn can properly account for tied edges", {
  # metadata are equal
  expect_identical_metadata(bruv_no_ties, bruv_ties)
  
  # There will be four tied edges, but only 3 untied.
  expect_ecount(bruv_no_ties, 3L)
  expect_ecount(bruv_ties, 4L)
  
  # vertex attributes come with pie
  expect_vertex_attr(bruv_no_ties, pienames)
  expect_vertex_attr(bruv_ties, pienames)
  
  # both graphs will have the same number of vertices
  expect_vcount(bruv_no_ties, 4L)
  expect_vcount(bruv_ties, 4L)
  
  # Adding a loop will shrink the distance between nodes 1 and 4
  expect_distance_table(bruv_no_ties, 3:1)
  expect_distance_table(bruv_ties, c(4, 2))
})

test_that("poppr.msn can properly account for tied edges", {
  # metadata are equal
  expect_identical_metadata(pmsn_no_ties, pmsn_ties)
  
  # There will be four tied edges, but only 3 untied.
  expect_ecount(pmsn_no_ties, 3L)
  expect_ecount(pmsn_ties, 4L)
  
  # vertex attributes come with pie
  expect_vertex_attr(pmsn_no_ties, pienames)
  expect_vertex_attr(pmsn_ties, pienames)
  
  # both graphs will have the same number of vertices
  expect_vcount(pmsn_no_ties, 4L)
  expect_vcount(pmsn_ties, 4L)
  
  # Adding a loop will shrink the distance between nodes 1 and 4
  expect_distance_table(pmsn_no_ties, 3:1)
  expect_distance_table(pmsn_ties, c(4, 2))
})

test_that("bruvo.msn can work with single populations", {
  # metadata are equal
  expect_identical_metadata(sbruv_no_ties, sbruv_ties)
  
  # There will be four tied edges, but only 3 untied.
  expect_ecount(sbruv_no_ties, 3L)
  expect_ecount(sbruv_ties, 4L)
  
  # both graphs will have the same number of vertices
  expect_vcount(sbruv_no_ties, 4L)
  expect_vcount(sbruv_ties, 4L)
  
  # vertex attributes don't come with pie :(
  expect_vertex_attr(sbruv_no_ties, nopienames)
  expect_vertex_attr(sbruv_ties, nopienames)
  
  # Adding a loop will shrink the distance between nodes 1 and 4
  expect_distance_table(sbruv_no_ties, 3:1)
  expect_distance_table(sbruv_ties, c(4, 2))
})

test_that("poppr.msn can work with single populations", {
  # metadata are equal
  expect_identical_metadata(spmsn_no_ties, spmsn_ties)
  
  # There will be four tied edges, but only 3 untied.
  expect_ecount(spmsn_no_ties, 3L)
  expect_ecount(spmsn_ties, 4L)
  
  # vertex attributes don't come with pie :(
  expect_vertex_attr(spmsn_no_ties, nopienames)
  expect_vertex_attr(spmsn_ties, nopienames)
  
  # both graphs will have the same number of vertices
  expect_vcount(spmsn_no_ties, 4L)
  expect_vcount(spmsn_ties, 4L)
  
  # Adding a loop will shrink the distance between nodes 1 and 4
  expect_distance_table(spmsn_no_ties, 3:1)
  expect_distance_table(spmsn_ties, c(4, 2))
})


context("MSN and collapsed MLG tests")

gmsnt <- bruvo.msn(gend, replen = c(1, 1), threshold = 0.15)

test_that("Minimum spanning networks also collapse MLGs", {
  skip_on_cran()
  gend        <- as.genclone(gend)
  gend_single <- as.genclone(gend_single)
  
  # Adding the filter for testing
  mlg.filter(gend, dist = bruvo.dist, replen = c(1, 1)) <- 0.15
  mll(gend) <- "original"
  
  expect_vcount(gmsnt, 2)

  pgmsnt <- poppr.msn(gend, distmat = gend_bruvo, threshold = 0.15)
  mll(gend) <- "contracted"
  gmsnot  <- bruvo.msn(gend, replen = c(1, 1)) # no threshold supplied
  gmsnone <- bruvo.msn(gend, replen = c(1, 1), threshold = 0.3)
  expect_vcount(gmsnone, 1)
  gmsnall  <- bruvo.msn(gend, replen = c(1, 1), threshold = 0)
  expect_vcount(gmsnall, 4)

  expect_identical_vertex_attr(gmsnt, pgmsnt, "pie")
  expect_identical_vertex_attr(gmsnot, pgmsnt, "pie")
  
  expect_identical_vertex_attr(gmsnt, pgmsnt, "name")
  expect_identical_vertex_attr(gmsnot, pgmsnt, "name")
  
  expect_identical_edge_attr(gmsnt, pgmsnt, "weight")
  expect_identical_edge_attr(gmsnot, pgmsnt, "weight")
  
  mll(gend) <- "original"
  gmsn   <- bruvo.msn(gend, replen = c(1, 1), showplot = FALSE)
  expect_vcount(gmsn, 4)
  pgmsn  <- poppr.msn(gend, distmat = gend_bruvo, showplot = FALSE)

  expect_identical_vertex_attr(gmsn, pgmsn, "pie")
  expect_identical_vertex_attr(gmsn, gmsnall, "pie")
  
  expect_identical_vertex_attr(gmsn, pgmsn, "name")
  expect_identical_vertex_attr(gmsn, gmsnall, "name")
  
  expect_identical_edge_attr(gmsn, pgmsn, "weight")
  expect_identical_edge_attr(gmsn, gmsnall, "weight")
  
})

test_that("Minimum spanning networks can collapse MLGs with single populations", {
  skip_on_cran()
  sgmsnt  <- bruvo.msn(gend_single, replen = c(1, 1), threshold = 0.15)
  psgmsnt <- poppr.msn(gend_single, distmat = gend_bruvo, threshold = 0.15)

  expect_identical_vertex_attr(sgmsnt, psgmsnt, "pie")
  expect_identical_vertex_attr(sgmsnt, psgmsnt, "name")
  expect_identical_edge_attr(sgmsnt, psgmsnt, "weight")
  
  sgmsn   <- bruvo.msn(gend_single, replen = c(1, 1), showplot = FALSE)
  psgmsn  <- poppr.msn(gend_single, distmat = gend_bruvo, showplot = FALSE)

  expect_identical_vertex_attr(sgmsn, psgmsn, "pie")
  expect_identical_vertex_attr(sgmsn, psgmsn, "name")
  expect_identical_edge_attr(sgmsn, psgmsn, "weight")

  expect_vcount(sgmsnt, 2)
  expect_vcount(sgmsn, 4)

  expect_output(plot_poppr_msn(gend, gmsnt, palette = "cm.colors"), NA)
  expect_output(plot_poppr_msn(gend_single, sgmsnt, palette = "cm.colors"), NA)
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

context("minimum spanning network subset populations")

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

context("custom MLLs and minimum spanning networks")

mll.custom(pc) <- LETTERS[mll(pc)]
mll.levels(pc)[mll.levels(pc) == "Q"] <- "M"
mll(pc) <- "custom"

test_that("msn works with custom MLLs", {
  skip_on_cran()
  expect_error(pcmsn <- bruvo.msn(pc, replen = rep(1, 10)), NA)
  expect_equivalent(sort(unique(igraph::V(pcmsn$graph)$label)), sort(mll.levels(pc)))
  expect_error(plot_poppr_msn(pc, pcmsn), NA)
  expect_error(plot_poppr_msn(pc, pcmsn, mlg = TRUE), NA)
})

context("Minimum spanning network aesthetics")

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
