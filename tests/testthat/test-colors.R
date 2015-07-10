context("Graph color conversion tests")

test_that("graph colors are correctly converted", {
  data(partial_clone)
  pc <- partial_clone[1:11]
  g <- bruvo.msn(pc, replen = rep(1, 10), showplot = FALSE)
  gu <- poppr:::update_poppr_graph(g, "rainbow")
  gpc <- igraph::V(g$graph)$pie.color
  gupc <- igraph::V(gu$graph)$pie.color
  expect_equal(names(unlist(gpc)), names(unlist(gupc)))
  expect_true(all(unique(unlist(gupc)) %in% rainbow(4)))
  expect_equivalent(gupc[[4]], rainbow(4)[3:4])
})
