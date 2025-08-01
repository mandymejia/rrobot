#' Compute SI Boot Thresholds for Outlier Detection
#'
#' Computes a robust distance (RD)–based threshold using single-imputed data
#' followed by bootstrap resampling over clean (included) indices. Returns the confidence interval
#' bounds of the bootstrapped 99th percentiles.
#'
#' @param RD_org_obj A list from `comp_RD()` on the original data.
#' @param imp_data A numeric matrix (T × Q) of single-imputed data.
#' @param B Integer; number of bootstrap resamples (default = 1000).
#' @param alpha Significance level for quantile thresholding (default = 0.01).
#' @param boot_quant Central quantile range for CI (default = 0.95, for 95% CI).
#' @param verbose Logical; whether to print progress.
#'
#' @return A list with:
#' \describe{
#'   \item{quant99}{Vector of 99th quantiles of RD for each bootstrap sample.}
#'   \item{LB_CI}{Lower bound of the confidence interval for the 99th quantiles.}
#'   \item{UB_CI}{Upper bound of the confidence interval for the 99th quantiles.}
#' }
#' @export
SI_boot <- function(RD_org_obj, imp_data,
                    B = 1000, alpha = 0.01, boot_quant = 0.95,
                    verbose = FALSE) {
  cov_mcd <- RD_org_obj$S_star
  ind_incld <- RD_org_obj$ind_incld
  n_time <- nrow(imp_data)
  
  quant99 <- numeric(B)
  cutoff_q <- 1 - alpha
  
  if (verbose) message("Running ", B, " bootstrap resamples...")
  
  for (b in seq_len(B)) {
    boot_idx <- sample(ind_incld, size = length(ind_incld), replace = TRUE)
    
    mu_boot <- colMeans(imp_data[boot_idx, , drop = FALSE])
    mu_mat <- matrix(mu_boot, nrow = n_time, ncol = ncol(imp_data), byrow = TRUE)
    
    RD_boot <- rowSums(((imp_data - mu_mat) %*% expm::sqrtm(solve(cov_mcd)))^2)
    quant99[b] <- quantile(RD_boot, cutoff_q, na.rm = TRUE)
    
    if (verbose && b %% 10 == 0) message("Bootstrap ", b, "/", B, " complete.")
  }
  
  lower_p <- (1 - boot_quant) / 2
  upper_p <- 1 - lower_p
  LB_CI <- quantile(quant99, probs = lower_p, na.rm = TRUE)
  UB_CI <- quantile(quant99, probs = upper_p, na.rm = TRUE)
  
  if (verbose) {
    message(round(boot_quant * 200), "% CI of 99th quantiles: [",
            round(LB_CI, 3), ", ", round(UB_CI, 3), "]")
  }
  
  return(list(
    quant99 = quant99,
    LB_CI = LB_CI,
    UB_CI = UB_CI
  ))
}
