#' Compute SI Boot Thresholds for Outlier Detection
#'
#' Computes a robust distance (RD)–based threshold using single-imputed data
#' followed by bootstrap resampling. Returns the confidence interval bounds of
#' the bootstrapped 99th percentiles.
#'
#' @param imp_data A numeric matrix (n_time × Q) of imputed data (from Step 2).
#' @param cov_mcd Robust covariance matrix (Q × Q) from original data (e.g., RD_obj$S_star).
#' @param ind_incld Row indices used to recompute robust mean (e.g., RD_obj$ind_incld).
#' @param B Integer; number of bootstrap resamples. Default is 1000.
#' @param alpha Significance level for thresholding (e.g., 0.01 for 99th percentile).
#' @param boot_quant Central quantile range for CI (default = 0.95, for 95% CI).
#' @param verbose Logical; whether to print progress messages.
#'
#' @return A list with:
#' \describe{
#'   \item{quant99}{Vector of 99th quantiles of RD for each bootstrap sample.}
#'   \item{LB_CI}{Lower bound of the confidence interval for the 99th quantiles.}
#'   \item{UB_CI}{Upper bound of the confidence interval for the 99th quantiles.}
#' }
#' @export
SI_boot <- function(imp_data, cov_mcd, ind_incld,
                    B = 1000, alpha = 0.01, boot_quant = 0.95,
                    verbose = FALSE) {
  n_time <- nrow(imp_data)
  
  # Store 99th quantiles for each bootstrap sample
  quant99 <- numeric(B)
  
  if (verbose) message("Running ", B, " bootstrap resamples...")
  
  for (b in seq_len(B)) {
    # Bootstrap sample (with replacement)
    boot_idx <- sample(seq_len(n_time), size = n_time, replace = TRUE)
    boot_sample <- imp_data[boot_idx, , drop = FALSE]
    
    # Compute RD using fixed covariance and robust mean from boot sample
    RD_obj <- comp_RD(
      data_matrix = boot_sample,
      mode = "manual",
      cov_mcd = cov_mcd,
      ind_incld = ind_incld
    )
    
    # Store 99th percentile
    quant99[b] <- quantile(RD_obj$RD, 1 - alpha, na.rm = TRUE)
    
    if (verbose && b %% 10 == 0) message("Bootstrap ", b, "/", B, " complete.")
  }
  
  # Compute CI bounds
  lower_p <- (1 - boot_quant) / 2
  upper_p <- 1 - lower_p
  LB_CI <- quantile(quant99, probs = lower_p, na.rm = TRUE)
  UB_CI <- quantile(quant99, probs = upper_p, na.rm = TRUE)
  
  if (verbose) {
    message(round(boot_quant * 100), "% CI of 99th quantiles: [",
            round(LB_CI, 3), ", ", round(UB_CI, 3), "]")
  }
  
  return(list(
    quant99 = quant99,
    LB_CI = LB_CI,
    UB_CI = UB_CI
  ))
}
