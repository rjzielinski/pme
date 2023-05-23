test_that("norm computation works", {
  expect_equal(norm_euclidean(1:10), norm(matrix(1:10, ncol = 1), type = "f"))
})
