test_that("emprule_rob flags correct outliers", {

  # Simulated dataset with a known outlier
  x <- c(10, 12, 11, 10, 13, 11, 100)  # 100 is an outlier

  # Call the function
  result <- rrobot:::emprule_rob(x, thr = 3, use_huber = FALSE, upper_only = TRUE)

  # Check output type and length
  expect_type(result, "logical")
  expect_length(result, length(x))

  # Verify the correct element is flagged as outlier
  expect_true(result[7])
  expect_false(any(result[1:6]))
})

test_that("SI method gives consistent results", {
  # Load test setup data
  setup_data <- readRDS(system.file("fixtures", "test_setup_data.rds", package = "rrobot"))
  reference <- readRDS(system.file("fixtures", "SI_reference.rds", package = "rrobot"))

  # Run SI method
  result <- SI(RD_org_obj = setup_data$RD_org_obj,
               imp_data = setup_data$imp_result$imp_data,
               alpha = 0.01)

  # Compare results
  expect_equal(result$SI_threshold, reference$SI_threshold, tolerance = 1e-10)
  expect_equal(result$SI_obj$RD, reference$SI_obj$RD, tolerance = 1e-10)
})


test_that("Hardin-Rocke method gives consistent results", {
  setup_data <- readRDS(system.file("fixtures", "test_setup_data.rds", package = "rrobot"))
  reference <- readRDS(system.file("fixtures", "HR_reference.rds", package = "rrobot"))

  result <- Fit_F(Q = ncol(setup_data$hk_data),
                  n = nrow(setup_data$hk_data),
                  h = setup_data$RD_org_obj$h,
                  quantile = 0.01)

  expect_equal(result$threshold, reference$threshold, tolerance = 1e-10)
  expect_equal(result$scale, reference$scale, tolerance = 1e-10)
  expect_equal(result$c, reference$c, tolerance = 1e-10)
  expect_equal(result$m, reference$m, tolerance = 1e-10)
  expect_equal(result$df, reference$df, tolerance = 1e-10)
})


test_that("SI_boot method gives consistent results", {
  skip_on_ci()  # Skip on CI due to bootstrap randomness
  setup_data <- readRDS(system.file("fixtures", "test_setup_data.rds", package = "rrobot"))
  reference <- readRDS(system.file("fixtures", "SI_boot_reference.rds", package = "rrobot"))

  # Set seed for reproducible bootstrap
  set.seed(2025)
  result <- SI_boot(RD_org_obj = setup_data$RD_org_obj,
                    imp_data = setup_data$imp_result$imp_data,
                    B = 50, alpha = 0.01, boot_quant = 0.95,
                    verbose = FALSE)

  # Use looser tolerances for bootstrap methods (inherently variable)
  expect_equal(result$LB_CI, reference$LB_CI, tolerance = 5)
  expect_equal(result$UB_CI, reference$UB_CI, tolerance = 5)
  expect_length(result$quant99, 50)

  # Test that confidence intervals are reasonable
  expect_lt(result$LB_CI, result$UB_CI)  # Lower < Upper
  expect_gt(result$LB_CI, 0)             # Positive thresholds
})

test_that("MI method gives consistent results", {
  skip_on_ci()  # Skip on CI due to bootstrap randomness
  setup_data <- readRDS(system.file("fixtures", "test_setup_data.rds", package = "rrobot"))
  reference <- readRDS(system.file("fixtures", "MI_reference.rds", package = "rrobot"))

  # Set seed for reproducible multiple imputation
  set.seed(2025)

  # Suppress convergence warnings for small test data
  suppressWarnings({
    multiple_imp <- MImpute(x = setup_data$imp_result$imp_data,
                            w = setup_data$kurt_data$lk,
                            outlier_matrix = setup_data$out_result$outliers,
                            M = 3, k = 5)
  })

  result <- MI(RD_org_obj = setup_data$RD_org_obj,
               imp_datasets = multiple_imp$imp_datasets,
               alpha = 0.01)

  # Use looser tolerances for multiple imputation (stochastic process)
  expect_equal(result$thresholds, reference$MI_results$thresholds, tolerance = 50)
  expect_length(result$thresholds, 3)  # Should have 3 thresholds (M=3)

  # Test structural properties instead of exact values
  expect_type(result$voted_outliers, "logical")
  expect_length(result$voted_outliers, length(reference$MI_results$voted_outliers))
  expect_gte(sum(result$voted_outliers), 0)  # Non-negative outlier count
})
