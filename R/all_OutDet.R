#' Apply All Outlier Detection Threshold Methods on One Dataset
#'
#' Performs ICA-based kurtosis split, univariate outlier detection and imputation,
#' then applies 5 outlier thresholding methods: HR, SI, SI_boot, MI, MI_boot.
#'
#' @param name Optional label for logging.
#' @param data_matrix T × V fMRI time series matrix.
#' @param alpha Quantile level for thresholding (default = 0.01).
#' @param B Number of bootstrap samples (default = 1000).
#' @param M Number of multiply imputed datasets (default = 50).
#' @param k Number of perturbation cycles per imputation (default = 100).
#'
#' @return A list of threshold values and intermediate diagnostics.
#' @export
all_OutDet <- function(name = NULL, data_matrix, alpha = 0.01, B = 1000, M = 50, k = 100) {
  if (!is.null(name)) message("Processing: ", name)
  
  # STEP 0: ICA-based Kurtosis Split
  kurt_data <- ICA_extract_kurt(time_series = data_matrix)
  x <- kurt_data$hk
  w <- kurt_data$lk
  Q <- ncol(x)
  
  # STEP 1: Univariate Temporal Outlier Detection
  out_result <- univOut(hk_data = x, cutoff = 3, trans = "SHASH")
  
  # STEP 2: Univariate Temporal Imputation
  imp_result <- impTemp_univOut(x = x, outlier_mask = out_result$outliers)
  x_imp <- imp_result$imp_data
  
  # STEP 3: Robust Distance on Original Data
  RD_org_obj <- comp_RD(data_matrix = x, mode = "auto")
  
  # STEP 4: Hard Rejection (HR) Threshold via F-approximation
  n <- length(RD_org_obj$RD)
  h <- length(RD_org_obj$ind_incld)
  fitF_obj <- Fit_F(Q = Q, n = n, h = h, quantile = alpha)
  q_HR <- fitF_obj$threshold * fitF_obj$scale
  
  # STEP 5: Multiple Imputation
  multiple_imp <- MImpute(
    x = x_imp,
    w = w,
    outlier_matrix = out_result$outliers,
    M = M,
    k = k,
    seed = 2025
  )
  
  # STEP 6: Apply 4 Robust Detection Methods
  SI_result      <- SI(RD_org_obj = RD_org_obj, imp_data = x_imp, alpha = alpha)
  SI_boot_result <- SI_boot(RD_org_obj = RD_org_obj, imp_data = x_imp, 
                            B = B, alpha = alpha, boot_quant = 0.95)
  MI_result      <- MI(RD_org_obj = RD_org_obj, imp_datasets = multiple_imp$imp_datasets, 
                       alpha = alpha)
  MI_boot_result <- MI_boot(RD_org_obj = RD_org_obj, imp_datasets = multiple_imp$imp_datasets, 
                            B = B, alpha = alpha, boot_quant = 0.95)
  
  # STEP 7: Return Summary
  return(list(
    HR        = q_HR,
    SI        = SI_result,
    SI_boot   = SI_boot_result,
    MI        = MI_result,
    MI_vote   = MI_result$voted_outliers,
    MI_boot   = MI_boot_result,
    name      = name,
    RD_org    = RD_org_obj
  ))
}
