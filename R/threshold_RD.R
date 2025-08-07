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
threshold_RD <- function(name = NULL, data_matrix, alpha = 0.01, B = 1000, M = 100, k = 200,
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
