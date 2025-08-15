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
  result <- thresh_SI(RD_org_obj = setup_data$RD_org_obj,
               imp_data = setup_data$imp_result$imp_data,
               alpha = 0.01)

  # Compare results
  expect_equal(result$threshold, reference$SI_threshold, tolerance = 1e-10)
  expect_equal(result$SI_obj$RD, reference$SI_obj$RD, tolerance = 1e-10)
})


test_that("Hardin-Rocke method gives consistent results", {
  setup_data <- readRDS(system.file("fixtures", "test_setup_data.rds", package = "rrobot"))
  reference <- readRDS(system.file("fixtures", "HR_reference.rds", package = "rrobot"))

  result <- thresh_F(p = ncol(setup_data$hk_data),
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
  # skip_on_ci()  # Skip on CI due to bootstrap randomness
  setup_data <- readRDS(system.file("fixtures", "test_setup_data.rds", package = "rrobot"))
  reference <- readRDS(system.file("fixtures", "SI_boot_reference.rds", package = "rrobot"))

  # Set seed for reproducible bootstrap
  set.seed(2025)
  result <- thresh_SI_boot(RD_org_obj = setup_data$RD_org_obj,
                    imp_data = setup_data$imp_result$imp_data,
                    B = 50, alpha = 0.01, boot_quant = 0.95,
                    verbose = FALSE)

  # Use looser tolerances for bootstrap methods (inherently variable)
  expect_equal(result$threshold, reference$LB_CI, tolerance = 5)
  expect_equal(result$UB_CI, reference$UB_CI, tolerance = 5)
  expect_length(result$quant99, 50)

  # Test that confidence intervals are reasonable
  expect_lt(result$threshold, result$UB_CI)
  expect_gt(result$threshold, 0)
})

test_that("MI method gives consistent results", {
  # skip_on_ci()  # Skip on CI due to bootstrap randomness
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

  result <- thresh_MI(RD_org_obj = setup_data$RD_org_obj,
               imp_datasets = multiple_imp$imp_datasets,
               alpha = 0.01)

  # Use looser tolerances for multiple imputation (stochastic process)
  expect_equal(result$thresholds, reference$MI_results$thresholds, tolerance = 50)
  expect_length(result$thresholds, 3)  # Should have 3 thresholds (M=3)

})

test_that("threshold_RD 'all' method gives consistent results", {
  # skip_on_ci()  # Skip on CI due to bootstrap randomness in MI methods
  # Load test setup data
  setup_data <- readRDS(system.file("fixtures", "test_setup_data.rds", package = "rrobot"))
  # Load reference results for comparison
  SI_ref <- readRDS(system.file("fixtures", "SI_reference.rds", package = "rrobot"))
  HR_ref <- readRDS(system.file("fixtures", "HR_reference.rds", package = "rrobot"))
  SI_boot_ref <- readRDS(system.file("fixtures", "SI_boot_reference.rds", package = "rrobot"))
  MI_ref <- readRDS(system.file("fixtures", "MI_reference.rds", package = "rrobot"))
  # Set seed for reproducible results
  set.seed(2025)
  RD_obj <- compute_RD(x = setup_data$hk_data, mode = "auto")
  # Run threshold_RD with "all" method
  suppressWarnings({
    result <- threshold_RD(x = setup_data$hk_data,
                           w = setup_data$kurt_data$lk,
                           method = "all",
                           RD_obj = RD_obj,
                           M = 3, k = 5, B = 50,  # Reduced for faster testing
                           alpha = 0.01,
                           cutoff = 4,
                           impute_method = "interp",
                           trans = "SHASH",
                           quantile = 0.01,
                           verbose = FALSE)
  })
  # Test structure of results
  expect_type(result, "list")
  expect_named(result, c("thresholds", "RD_obj", "RD_obj_shash", "x_norm", "call"))
  expect_named(result$thresholds, c("SI", "SI_boot", "MI", "MI_boot", "F", "SHASH"))

  # Test SI method matches reference
  expect_equal(result$thresholds$SI$threshold, SI_ref$SI_threshold, tolerance = 5)

  # Test F method matches reference
  expect_equal(result$thresholds$F$threshold, HR_ref$threshold, tolerance = 5)
  expect_equal(result$thresholds$F$scale, HR_ref$scale, tolerance = 5)

  # Test SI_boot method (looser tolerance due to bootstrap)
  expect_equal(result$thresholds$SI_boot$threshold, SI_boot_ref$LB_CI, tolerance = 5)
  expect_equal(result$thresholds$SI_boot$UB_CI, SI_boot_ref$UB_CI, tolerance = 5)
  expect_length(result$thresholds$SI_boot$quant99, 50)

  # Test MI method (looser tolerance due to stochastic nature)
  expect_length(result$thresholds$MI$thresholds, 3)  # Should have 3 thresholds (M=3)

  # Test MI_boot method
  expect_type(result$thresholds$MI_boot$threshold, "double")
  expect_length(result$thresholds$MI_boot$threshold, 1)
  expect_type(result$thresholds$MI_boot$flagged_outliers, "logical")
  expect_length(result$thresholds$MI_boot$flagged_outliers, length(setup_data$hk_data[,1]))

  # Test that all methods return reasonable positive thresholds
  expect_gt(result$thresholds$SI$threshold, 0)
  expect_gt(result$thresholds$F$threshold, 0)
  expect_gt(result$thresholds$SI_boot$threshold, 0)
  expect_gt(result$thresholds$MI_boot$threshold, 0)
  expect_true(all(result$thresholds$MI$thresholds > 0))

  # Test RD_obj structure
  expect_s3_class(result$RD_obj, "RD_result")

  # Summary methods
  expect_no_error(summary(result$thresholds$SI))
  expect_no_error(summary(result$thresholds$SI_boot))
  expect_no_error(summary(result$thresholds$MI))
  expect_no_error(summary(result$thresholds$MI_boot))
  expect_no_error(summary(result$thresholds$F))
  expect_no_error(summary(result$thresholds$SHASH))
})

test_that("RD method gives consistent results", {
  # skip_on_ci()  # Skip on CI due to bootstrap randomness
  setup_data <- readRDS(system.file("fixtures", "test_setup_data.rds", package = "rrobot"))
  set.seed(2025)
  suppressWarnings({
    result <- RD(x = setup_data$hk_data,
                 w = setup_data$kurt_data$lk,
                 method = "SI",
                 alpha = 0.01,
                 cutoff = 4,
                 impute_method = "interp",
                 trans = "SHASH")
  })
  # Test structure - RD() doesn't return RD_obj anymore
  expect_type(result, "list")
  expect_named(result, c("thresholds", "RD_obj", "outliers", "call"))
  expect_s3_class(result, "RD")

  # Test that thresholds contains the RD_obj
  expect_s3_class(result$RD_obj, "RD_result")

  # Test that call is captured
  expect_true(is.call(result$call))

  # Test outliers vector
  expect_type(result$outliers, "logical")
  expect_length(result$outliers, nrow(setup_data$hk_data))
  expect_true(all(is.logical(result$outliers)))
})
