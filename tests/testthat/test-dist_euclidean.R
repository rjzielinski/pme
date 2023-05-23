test_that("distance computation works", {
  expect_equal(dist_euclidean(1:10, 91:100), norm(matrix(1:10 - 91:100, ncol = 1), type = "f"))
})
