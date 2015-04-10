context("Errors and Exceptions")

test_that("bruvo.msn throws errors with incorrect add and loss arguments", {
  data("partial_clone", package = "poppr")
  reps <- rep(1, 10)
  expect_error(bruvo.msn(partial_clone[1:10], replen = reps, add = "1", loss = TRUE))
  expect_error(bruvo.msn(partial_clone[1:10], replen = reps, add = TRUE, loss = "TRUE"))
  expect_error(bruvo.msn(partial_clone[1:10], replen = reps, add = matrix(TRUE, 10, 10), loss = TRUE))
})