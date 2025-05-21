test_that("emprule_rob flags correct outliers", {

  # Simulated dataset with a known outlier
  x <- c(10, 12, 11, 10, 13, 11, 100)  # 100 is an outlier

  # Call the function
  result <- rrobot::emprule_rob(x, thr = 3, use_huber = FALSE, upper_only = TRUE)

  # Check output type and length
  expect_type(result, "logical")
  expect_length(result, length(x))

  # Verify the correct element is flagged as outlier
  expect_true(result[7])
  expect_false(any(result[1:6]))
})
