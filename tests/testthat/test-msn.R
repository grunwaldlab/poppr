context("Minimum spanning network functions")

test_that("Minimum spanning networks can properly account for tied edges", {

set.seed(9005)

gend <-new("genind"
           , tab = structure(c(0.5, 0.5, 0, 0, 0.5, 0.5, 1, 1, 0, 0.5, 0.5, 0, 1, 
                               0.5, 0.5, 1), .Dim = c(4L, 4L), .Dimnames = list(c("1", "2", 
                                                                                  "3", "4"), c("L1.1", "L1.2", "L2.1", "L2.2")))
           , loc.names = structure(c("loc-1", "loc-2"), .Names = c("L1", "L2"))
           , loc.fac = structure(c(1L, 1L, 2L, 2L), .Label = c("L1", "L2"), class = "factor")
           , loc.nall = structure(c(2L, 2L), .Names = c("L1", "L2"))
           , all.names = structure(list(L1 = structure(c("3", "4"), .Names = c("1", "2"
           )), L2 = structure(c("3", "4"), .Names = c("1", "2"))), .Names = c("L1", 
                                                                              "L2"))
           , call = NULL
           , ind.names = structure(c("", "", "", ""), .Names = c("1", "2", "3", "4"))
           , pop = structure(1:4, .Label = c("P1", "P2", "P3", "P4"), class = "factor")
           , pop.names = structure(c("1", "2", "3", "4"), .Names = c("P1", "P2", "P3", 
                                                                     "P4"))
           , ploidy = 2L
           , type = "codom"
           , other = NULL
) 

set.seed(9005)

gend_single <- new("genind"
    , tab = structure(c(0.5, 0.5, 0, 0, 0.5, 0.5, 1, 1, 0, 0.5, 0.5, 0, 1, 
0.5, 0.5, 1), .Dim = c(4L, 4L), .Dimnames = list(c("1", "2", 
"3", "4"), c("L1.1", "L1.2", "L2.1", "L2.2")))
    , loc.names = structure(c("loc-1", "loc-2"), .Names = c("L1", "L2"))
    , loc.fac = structure(c(1L, 1L, 2L, 2L), .Label = c("L1", "L2"), class = "factor")
    , loc.nall = structure(c(2L, 2L), .Names = c("L1", "L2"))
    , all.names = structure(list(L1 = structure(c("3", "4"), .Names = c("1", "2"
)), L2 = structure(c("3", "4"), .Names = c("1", "2"))), .Names = c("L1", 
"L2"))
    , call = NULL
    , ind.names = structure(c("", "", "", ""), .Names = c("1", "2", "3", "4"))
    , pop = NULL
    , pop.names = NULL
    , ploidy = 2L
    , type = "codom"
    , other = NULL
)

no_ties <- structure(list(graph = structure(list(4, FALSE, c(3, 2, 3), c(0, 
1, 2), c(1, 0, 2), c(0, 1, 2), c(0, 0, 0, 1, 3), c(0, 1, 2, 3, 
3), list(c(1, 0, 1), structure(list(), .Names = character(0)), 
    structure(list(name = c("P1", "P2", "P3", "P4"), size = c(1L, 
    1L, 1L, 1L), shape = c("pie", "pie", "pie", "pie"), pie = list(
        structure(1L, .Names = "1"), structure(1L, .Names = "2"), 
        structure(1L, .Names = "3"), structure(1L, .Names = "4")), 
        pie.color = list("#4C00FFFF", "#00E5FFFF", "#00FF4DFF", 
            "#FFFF00FF"), label = c("MLG.3", "MLG.4", "MLG.2", 
        "MLG.1")), .Names = c("name", "size", "shape", "pie", 
    "pie.color", "label")), structure(list(weight = c(0.125, 
    0.125, 0.125), color = c("#434343", "#434343", "#434343"), 
        width = c(8, 8, 8)), .Names = c("weight", "color", "width"
    )))), class = "igraph"), populations = structure(c("1", "2", 
"3", "4"), .Names = c("P1", "P2", "P3", "P4")), colors = c("#4C00FFFF", 
"#00E5FFFF", "#00FF4DFF", "#FFFF00FF")), .Names = c("graph", 
"populations", "colors"))

no_ties_single <- structure(list(graph = structure(list(4, FALSE, c(3, 2, 3), c(0, 
1, 2), c(1, 0, 2), c(0, 1, 2), c(0, 0, 0, 1, 3), c(0, 1, 2, 3, 
3), list(c(1, 0, 1), structure(list(), .Names = character(0)), 
    structure(list(size = c(1L, 1L, 1L, 1L), color = c("#4C00FFFF", 
    "#4C00FFFF", "#4C00FFFF", "#4C00FFFF"), label = c("MLG.3", 
    "MLG.4", "MLG.2", "MLG.1")), .Names = c("size", "color", 
    "label")), structure(list(weight = c(0.125, 0.125, 0.125), 
        color = c("#434343", "#434343", "#434343"), width = c(8, 
        8, 8)), .Names = c("weight", "color", "width")))), class = "igraph"), 
    populations = NA, colors = "#4C00FFFF"), .Names = c("graph", 
"populations", "colors"))

ties <- structure(list(graph = structure(list(4, FALSE, c(3, 2, 3, 1), 
    c(0, 1, 2, 0), c(3, 1, 0, 2), c(3, 0, 1, 2), c(0, 0, 1, 2, 
    4), c(0, 2, 3, 4, 4), list(c(1, 0, 1), structure(list(), .Names = character(0)), 
        structure(list(name = c("P1", "P2", "P3", "P4"), size = c(1L, 
        1L, 1L, 1L), shape = c("pie", "pie", "pie", "pie"), pie = list(
            structure(1L, .Names = "1"), structure(1L, .Names = "2"), 
            structure(1L, .Names = "3"), structure(1L, .Names = "4")), 
            pie.color = list("#4C00FFFF", "#00E5FFFF", "#00FF4DFF", 
                "#FFFF00FF"), label = c("MLG.3", "MLG.4", "MLG.2", 
            "MLG.1")), .Names = c("name", "size", "shape", "pie", 
        "pie.color", "label")), structure(list(weight = c(0.125, 
        0.125, 0.125, 0.125), color = c("#434343", "#434343", 
        "#434343", "#434343"), width = c(8, 8, 8, 8)), .Names = c("weight", 
        "color", "width")))), class = "igraph"), populations = structure(c("1", 
"2", "3", "4"), .Names = c("P1", "P2", "P3", "P4")), colors = c("#4C00FFFF", 
"#00E5FFFF", "#00FF4DFF", "#FFFF00FF")), .Names = c("graph", 
"populations", "colors"))

ties_single <- structure(list(graph = structure(list(4, FALSE, c(3, 2, 3), c(0, 
1, 2), c(1, 0, 2), c(0, 1, 2), c(0, 0, 0, 1, 3), c(0, 1, 2, 3, 
3), list(c(1, 0, 1), structure(list(), .Names = character(0)), 
    structure(list(size = c(1L, 1L, 1L, 1L), color = c("#4C00FFFF", 
    "#4C00FFFF", "#4C00FFFF", "#4C00FFFF"), label = c("MLG.3", 
    "MLG.4", "MLG.2", "MLG.1")), .Names = c("size", "color", 
    "label")), structure(list(weight = c(0.125, 0.125, 0.125), 
        color = c("#434343", "#434343", "#434343"), width = c(8, 
        8, 8)), .Names = c("weight", "color", "width")))), class = "igraph"), 
    populations = NA, colors = "#4C00FFFF"), .Names = c("graph", 
"populations", "colors"))


# Test Bruvo.msn
set.seed(9005)
expect_that(bruvo.msn(gend, replen=c(1,1)), equals(no_ties))
set.seed(9005)
expect_that(bruvo.msn(gend, replen=c(1,1), include.ties = TRUE), equals(ties))

# Test poppr.msn
set.seed(9005)
expect_that(poppr.msn(gend, distmat=bruvo.dist(gend,replen=c(1,1))), equals(no_ties))
set.seed(9005)
expect_that(poppr.msn(gend, distmat=bruvo.dist(gend,replen=c(1,1)), include.ties = TRUE), equals(ties))

# Test both for single populations sets
set.seed(9005)
expect_that(bruvo.msn(gend_single, replen=c(1,1)), equals(no_ties_single))
set.seed(9005)
expect_that(bruvo.msn(gend_single, replen=c(1,1), include.ties = TRUE), equals(ties_single))
set.seed(9005)
expect_that(poppr.msn(gend_single, distmat=bruvo.dist(gend,replen=c(1,1))), equals(no_ties_single))
set.seed(9005)
expect_that(poppr.msn(gend_single, distmat=bruvo.dist(gend,replen=c(1,1)), include.ties = TRUE), equals(ties_single))



})
