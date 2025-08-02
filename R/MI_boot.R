#' Bootstrap-Based Outlier Detection via Multiple Imputation (MI_boot)
#'
#' Extends single imputation bootstrapping by using multiple imputation. For each of the M imputed datasets:
#' - Generate B bootstrap samples (with replacement) from the included (non-outlier) indices
#' - Compute robust distances (RD) for each sample
#' - Extract the (1 - alpha) quantile of RD for each sample
#'
#' This yields M × B threshold candidates. The lower bound of their (1 - boot_quant) confidence interval
#' is used as the final threshold. This is applied to the RD of the original data.
#'
#' @param RD_org_obj A list output from `comp_RD()` on the original hk_data.
#' @param imp_datasets A list of M numeric matrices (T × Q); multiply imputed datasets.
#' @param B Number of bootstrap samples per imputed dataset (default = 1000).
#' @param alpha Significance level for quantile thresholding (default = 0.01).
#' @param boot_quant Confidence level for CI (default = 0.95).
#' @param verbose Logical; if TRUE, prints progress.
#'
#' @return A list with:
#' \describe{
#'   \item{thresholds_all}{Vector of M×B thresholds from each bootstrap sample.}
#'   \item{final_threshold}{Lower bound of CI across thresholds.}
#'   \item{flagged_outliers}{Logical vector of outliers based on final threshold.}
#' }
#' @export
MI_boot <- function(RD_org_obj, imp_datasets, B = 1000, alpha = 0.01, boot_quant = 0.95, verbose = FALSE) {
  M <- length(imp_datasets)
  n_time <- length(RD_org_obj$RD)
  Q <- ncol(imp_datasets[[1]])
  thresholds_all <- numeric(M * B)

  cov_mcd <- RD_org_obj$S_star
  ind_incld <- RD_org_obj$ind_incld
  RD_orig <- RD_org_obj$RD
  cutoff_q <- 1 - alpha

  invcov_sqrt <- RD_org_obj$invcov_sqrt

  idx <- 1
  for (m in seq_len(M)) {
    imp_data <- imp_datasets[[m]]
    stopifnot(nrow(imp_data) == n_time, ncol(imp_data) == Q)

    for (b in seq_len(B)) {
      boot_indices <- sample(ind_incld, size = length(ind_incld), replace = TRUE)

      mu_boot <- colMeans(imp_data[boot_indices, , drop = FALSE])
      mu_mat <- matrix(mu_boot, nrow = n_time, ncol = Q, byrow = TRUE)

      RD_boot <- rowSums(((imp_data - mu_mat) %*% invcov_sqrt)^2)
      thresholds_all[idx] <- quantile(RD_boot, cutoff_q, na.rm = TRUE)
      idx <- idx + 1

      if (verbose && b %% 100 == 0) message(sprintf("MI %d, Bootstrap %d", m, b))
    }
    if (verbose) message(sprintf("Completed MI %d of %d", m, M))
  }

  lb_ci <- quantile(thresholds_all, probs = 1 - boot_quant, na.rm = TRUE)
  flagged_outliers <- RD_orig > lb_ci

  result <- list(
    thresholds_all = thresholds_all,
    final_threshold = unname(lb_ci),
    flagged_outliers = flagged_outliers
  )

  class(result) <- "MI_boot_result"
  return(result)
}
