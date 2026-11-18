test_that("emprule_rob flags correct outliers", {
  x <- c(10, 12, 11, 10, 13, 11, 100)
  result <- rrobot:::emprule_rob(x, thr = 3, tail = "upper")

  expect_type(result, "logical")
  expect_length(result, length(x))
  expect_true(result[7])
  expect_false(any(result[1:6]))
})

test_that("univOut records the call", {
  result <- rrobot:::univOut(matrix(rnorm(20), ncol = 1), method = "robZ")

  expect_true(is.call(result$call))
  expect_identical(result$call$method, "robZ")
})

test_that("SHASH transforms round trip", {
  x <- rnorm(50)
  z <- normal_to_SHASH(x, mu = 1, sigma = 2, nu = 0.2, tau = 1.1)
  recovered <- SHASH_to_normal(z, mu = 1, sigma = 2, nu = 0.2, tau = 1.1)

  expect_equal(recovered, x, tolerance = 1e-10)
})
