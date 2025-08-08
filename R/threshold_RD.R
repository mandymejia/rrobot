#' Apply Outlier Detection Thresholding Methods
#'
#' Applies five robust thresholding methods (HR, SI, SI_boot, MI, MI_boot)
#' after ICA-based univariate outlier detection and imputation.
#'
#' @param name Optional label for logging.
#' @param data_matrix T × V fMRI time series matrix.
#' @param alpha Quantile level for thresholding (default = 0.01).
#' @param B Number of bootstrap samples (default = 1000).
#' @param M Number of multiply imputed datasets (default = 50).
#' @param k Number of perturbation cycles per imputation (default = 100).
#' @param trans Character; transformation method for outlier detection (default = "SHASH").
#' @param cutoff Numeric; MAD-based outlier cutoff (default = 3).
#'
#' @return A list containing threshold values and intermediate diagnostics.
#' @export
threshold_RD_old <- function(name = NULL, data_matrix, alpha = 0.01, B = 1000, M = 100, k = 200,
                       trans = "SHASH", cutoff = 3) {
  if (!is.null(name)) message("Processing: ", name)
  # STEP 0: ICA-based Kurtosis Split
  # kurt_data <- ICA_extract_kurt(time_series = data_matrix)
  # x <- kurt_data$hk
  # w <- kurt_data$lk
  x <- data_matrix
  w <- data_matrix
  Q <- ncol(x)
  # STEP 1: Univariate Temporal Outlier Detection
  out_result_x <- univOut(hk_data = x, cutoff = cutoff, trans = trans)
  imp_result_x <- impTemp_univOut(x = x, outlier_mask = out_result_x$outliers)
  x_imp <- imp_result_x$imp_data
  # STEP 3: Robust Distance on Original Data
  RD_org_obj <- RD(data_matrix = x, mode = "auto")
  # STEP 4: Hard Rejection (HR) Threshold via F-approximation
  n <- length(RD_org_obj$RD)
  h <- length(RD_org_obj$ind_incld)
  fitF_obj <- thresh_F(Q = Q, n = n, h = h, quantile = alpha)
  q_HR <- fitF_obj$threshold * fitF_obj$scale
  # STEP 5: Multiple Imputation
  # Note: w = NULL disables lk predictors
  multiple_imp <- MImpute(
    x = x_imp,
    w = w,
    outlier_matrix = out_result_x$outliers,
    M = M,
    k = k
  )
  # STEP 6: Apply 4 Robust Detection Methods
  SI_result      <- thresh_SI(RD_org_obj = RD_org_obj, imp_data = x_imp, alpha = alpha)
  SI_boot_result <- thresh_SI_boot(RD_org_obj = RD_org_obj, imp_data = x_imp,
                            B = B, alpha = alpha, boot_quant = 0.95)
  MI_result      <- thresh_MI(RD_org_obj = RD_org_obj, imp_datasets = multiple_imp$imp_datasets,
                       alpha = alpha)
  MI_boot_result <- thresh_MI_boot(RD_org_obj = RD_org_obj, imp_datasets = multiple_imp$imp_datasets,
                            B = B, alpha = alpha, boot_quant = 0.95)
  # STEP 7: Return Summary
  return(list(
    # Threshold values from all methods
    HR        = q_HR,
    SI        = SI_result,
    SI_boot   = SI_boot_result,
    MI        = MI_result,
    MI_boot   = MI_boot_result,
    # Optional: voting output from MI
    MI_vote   = MI_result$voted_outliers,
    # Metadata and diagnostics
    name      = name,
    RD_org    = RD_org_obj
  ))
}

#' Comprehensive Outlier Detection Using Robust Distance Thresholding
#'
#' Performs complete outlier detection workflow including univariate outlier detection,
#' imputation, robust distance computation, and threshold estimation using various methods.
#'
#' @param x A numeric matrix or data frame of dimensions T × Q (observations × variables).
#' @param w A numeric matrix (n_time × L) of low-kurtosis ICA components used as predictors.
#' @param threshold_method Character string specifying threshold method: "all", "SI", "SI_boot", "MI", "MI_boot", or "F".
#' @param RD_obj Optional pre-computed RD object from RD() function. If NULL, will compute RD internally.
#' @param mode Character string; either "auto" (default) to compute MCD internally or "manual" to use user-supplied values.
#' @param cov_mcd Optional covariance matrix (Q × Q); required in "manual" mode.
#' @param ind_incld Optional vector of row indices used to compute the robust mean; required in "manual" mode.
#' @param dist Logical; if TRUE, compute squared robust Mahalanobis distances for all observations.
#' @param cutoff A numeric value indicating how many MADs away from the median to flag as outliers. The default value is set to be 4.
#' @param trans Character string. One of "SHASH", "robust-YJ", or "robMAD".
#' @param M Number of multiply imputed datasets (default = 5).
#' @param k Number of perturbation cycles per imputation (default = 10).
#' @param alpha Significance level for threshold computation (default = 0.01).
#' @param quantile Numeric in (0,1) specifying the upper quantile for thresholding.
#' @param verbose Logical; if TRUE, print progress messages.
#' @param boot_quant Confidence level for CI (default = 0.95).
#' @param B Number of bootstrap samples per imputation (default = 1000).
#'
#' @export
threshold_RD <- function(x, w = NULL, threshold_method = c("all", "SI", "SI_boot", "MI", "MI_boot", "F"), RD_obj = NULL,
                                     # RD parameters
                                     mode = "auto", cov_mcd = NULL, ind_incld = NULL, dist = TRUE,
                                     # univOut parameters
                                     cutoff = 4,
                                     trans = "SHASH",
                                     # MImpute parameters
                                     M = 50,
                                     k = 100,
                                     # Threshold parameters
                                     alpha = 0.01,
                                     quantile = 0.01,
                                     verbose = FALSE,
                                     boot_quant = 0.95,
                                     B = 1000) {
  threshold_method <- match.arg(threshold_method)

  # Data pre-processing
  RD_obj <- if (is.null(RD_obj)) RD(data_matrix = x, mode = mode, cov_mcd = cov_mcd, ind_incld = ind_incld, dist = dist) else RD_obj

  out_result <- univOut(hk_data = x, cutoff = cutoff, trans = trans)
  imp_result <- impTemp_univOut(x = x, outlier_mask = out_result$outliers)


  if (threshold_method %in% c("all", "MI", "MI_boot")) {
    multiple_imp <- MImpute(x = imp_result$imp_data, w = w, outlier_matrix = out_result$outliers, M = M, k = k)
  }

  result <- switch(threshold_method,
                   "SI" = thresh_SI(RD_org_obj = RD_obj, imp_data = imp_result$imp_data, alpha = alpha),

                   "SI_boot" = thresh_SI_boot(RD_org_obj = RD_obj, imp_data = imp_result$imp_data,
                                              B = B, alpha = alpha, boot_quant = boot_quant, verbose = verbose),

                   "MI" = thresh_MI(RD_org_obj = RD_obj, imp_datasets = multiple_imp$imp_datasets, alpha = alpha),

                   "MI_boot" = thresh_MI_boot(RD_org_obj = RD_obj, imp_datasets = multiple_imp$imp_datasets,
                                              B = B, alpha = alpha, boot_quant = boot_quant, verbose = verbose),

                   "F" = thresh_F(Q = ncol(x), n = nrow(x), h = RD_obj$h, quantile = quantile),

                   "all" = list(
                     SI = thresh_SI(RD_org_obj = RD_obj, imp_data = imp_result$imp_data, alpha = alpha),
                     SI_boot = thresh_SI_boot(RD_org_obj = RD_obj, imp_data = imp_result$imp_data,
                                              B = B, alpha = alpha, boot_quant = boot_quant, verbose = verbose),
                     MI = thresh_MI(RD_org_obj = RD_obj, imp_datasets = multiple_imp$imp_datasets, alpha = alpha),
                     MI_boot = thresh_MI_boot(RD_org_obj = RD_obj, imp_datasets = multiple_imp$imp_datasets,
                                              B = B, alpha = alpha, boot_quant = boot_quant, verbose = verbose),
                     F = thresh_F(Q = ncol(x), n = nrow(x), h = RD_obj$h, quantile = quantile)
                   )
  )


  result
}
