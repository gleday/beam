context("mycontext")

data("TCPAprad")

test_that("bla", {
  fit <- beam(X = TCPAprad, type="both", verbose=FALSE)
  expect_equal(fit@dimX, c(164, 189))
})