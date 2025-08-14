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
#' @inheritParams RD_org_obj
#' @inheritParams imp_datasets
#' @inheritParams B
#' @inheritParams alpha
#' @inheritParams boot_quant
#' @inheritParams verbose
#'
#' @return A list with:
#' \describe{
#'   \item{thresholds_all}{Vector of M×B thresholds from each bootstrap sample.}
#'   \item{threshold}{Lower bound of CI across thresholds.}
#'   \item{flagged_outliers}{Logical vector of outliers based on final threshold.}
#'   \item{call}{The matched function call.}
#' }
#' @export
thresh_MI_boot <- function(RD_org_obj, imp_datasets, B = 1000, alpha = 0.01, boot_quant = 0.95, verbose = FALSE) {
  if (verbose) message("Running MI_boot method: ", length(imp_datasets)*B, " total bootstrap samples...")
  call <- match.call()

  M <- length(imp_datasets)
  n_time <- length(RD_org_obj$RD)
  p <- ncol(imp_datasets[[1]])
  thresholds_all <- numeric(M * B)

  RD_org <- RD_org_obj$RD
  ind_incld <- RD_org_obj$ind_incld
  ind_excld <- RD_org_obj$ind_excld
  invcov_sqrt <- RD_org_obj$invcov_sqrt
  cutoff_q <- 1 - alpha

  idx <- 1

  for (m in seq_len(M)) {
    imp_data <- imp_datasets[[m]]
    stopifnot(nrow(imp_data) == n_time, ncol(imp_data) == p)

    for (b in seq_len(B)) {
      # Resample both included and excluded indices
      boot_incld <- sample(ind_incld, size = length(ind_incld), replace = TRUE)
      boot_excld <- sample(ind_excld, size = length(ind_excld), replace = TRUE)

      # Combine and sort indices (sorting optional)
      boot_idx <- sort(c(boot_incld, boot_excld))

      # Bootstrap sample of the imputed dataset
      imp_boot <- imp_data[boot_idx, , drop = FALSE]

      # Compute mean and RD on the bootstrapped data
      # mu_boot <- colMeans(imp_boot)
      mu_boot <- RD_org_obj$xbar_star # Don't re-compute mean
      mu_mat <- matrix(mu_boot, nrow = nrow(imp_boot), ncol = p, byrow = TRUE)

      RD_boot <- rowSums(((imp_boot - mu_mat) %*% invcov_sqrt)^2)
      thresholds_all[idx] <- quantile(RD_boot, cutoff_q, na.rm = TRUE)
      idx <- idx + 1

      if (verbose && b %% 100 == 0) message(sprintf("MI %d, Bootstrap %d", m, b))
    }

    if (verbose) message(sprintf("Completed MI %d of %d", m, M))
  }


  lb_ci <- quantile(thresholds_all, probs = 1 - boot_quant, na.rm = TRUE)
  flagged_outliers <- RD_org > lb_ci

  result <- list(
    thresholds_all = thresholds_all,
    threshold = unname(lb_ci),
    flagged_outliers = flagged_outliers,
    call = call
  )

  class(result) <- "MI_boot_result"
  result

}
