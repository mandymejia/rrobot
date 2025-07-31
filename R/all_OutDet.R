#' Apply All Outlier Detection Threshold Methods on One Dataset
#'
#' @param name String. A label for the dataset (e.g., "abide1")
#' @param x High-kurtosis matrix (T × K)
#' @param w Low-kurtosis matrix (T × L)
#' @param Q Integer. Number of variables in x (usually ncol(x))
#' @param alpha Numeric. Quantile for threshold (default 0.01)
#'
#' @return A list containing thresholds for each method
all_OutDet <- function(name, x, w, Q, alpha = 0.01) {
  message("Processing: ", name)
  
  # Step 1: Compute RD for original data (needed for HR)
  RD_org_obj <- comp_RD(data_matrix = x, mode = "auto")
  
  # Step 2: HR threshold
  n <- length(RD_org_obj$RD)
  h <- length(RD_org_obj$ind_incld)
  fitF_obj <- Fit_F(Q = Q, n = n, h = h, quantile = alpha)
  q_HR <- fitF_obj$threshold * fitF_obj$scale
  
  # Step 3: Run all other imputation-based thresholds
  SI_result      <- SI(x = x, w = w, cutoff = 4, trans = "SHASH")
  SI_boot_result <- SI_boot(x = x, w = w, cutoff = 4, trans = "SHASH")
  MI_result      <- MI(x = x, w = w, M = 5, k = 10, cutoff = 4, trans = "SHASH")
  MI_boot_result <- MI_boot(x = x, w = w, B = 10, M = 5, k = 10, cutoff = 4, trans = "SHASH")
  
  list(
    HR       = q_HR,
    SI       = SI_result$threshold,
    SI_boot  = SI_boot_result$threshold,
    MI       = MI_result$thresholds,
    MI_vote  = MI_result$voted_cutoff_description,
    MI_boot  = MI_boot_result$threshold
  )
}
