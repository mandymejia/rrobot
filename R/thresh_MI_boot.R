#' Bootstrap-Based Outlier Detection via Multiple Imputation (MI_boot)
#'
#' Extends single imputation bootstrapping by using multiple imputation. For each of the M imputed datasets:
#' - Generate B bootstrap samples; included and excluded indices (from the robust subset)
#'   are resampled separately, with replacement
#' - Compute robust distances (RD) for each sample
#' - Extract the (1 - alpha) quantile of RD for each sample
#'
#' This yields M x B threshold candidates. The lower bound of their confidence interval
#' is used as the final threshold. This is applied to the RD of the original data.
#'
#' @inheritParams params_doc
#'
#' @return A list of class \code{"MI_boot_result"} with:
#' \describe{
#'   \item{thresholds}{Vector of M x B thresholds from each bootstrap sample.}
#'   \item{threshold}{Lower bound of CI across thresholds.}
#'   \item{flagged_outliers}{Integer vector of row indices from original data matrix that exceed the threshold.}
#'   \item{n_obs}{Number of observations in the original data.}
#'   \item{call}{The matched function call.}
#' }
#' @importFrom stats quantile mahalanobis
#' @keywords internal
thresh_MI_boot <- function(RD_org_obj, imp_datasets,
                           B = 1000, alpha = 0.01, boot_quant = 0.95,
                           verbose = FALSE) {
  call <- match.call()
  stopifnot(!is.null(RD_org_obj$RD), !is.null(RD_org_obj$S_star))
  if (verbose) message("Running MI_boot method: ", length(imp_datasets) * B, " total bootstrap samples...")

  M <- length(imp_datasets)
  n_time <- length(RD_org_obj$RD)
  p <- ncol(imp_datasets[[1]])
  thresholds_all <- numeric(M * B)

  RD_org    <- RD_org_obj$RD
  ind_incld <- RD_org_obj$ind_incld
  ind_excld <- RD_org_obj$ind_excld
  center0   <- RD_org_obj$xbar_star
  S0        <- RD_org_obj$S_star
  cutoff_q  <- 1 - alpha

  n_incld <- length(ind_incld)
  n_excld <- length(ind_excld)

  idx <- 1

  for (m in seq_len(M)) {
    imp_data <- imp_datasets[[m]]
    stopifnot(nrow(imp_data) == n_time, ncol(imp_data) == p)

    # Distances are row-wise: compute once per imputed dataset (original robust
    # center/covariance) and resample the vector, instead of recomputing for
    # every bootstrap sample. Gives identical values.
    RD_imp <- stats::mahalanobis(imp_data, center = center0, cov = S0)

    for (b in seq_len(B)) {
      # Resample both included and excluded indices
      # (sample.int avoids sample()'s 1:n behaviour when an index set has length 1)
      boot_idx <- c(ind_incld[sample.int(n_incld, replace = TRUE)],
                    ind_excld[sample.int(n_excld, replace = TRUE)])

      thresholds_all[idx] <- stats::quantile(RD_imp[boot_idx], cutoff_q, na.rm = TRUE, names = FALSE)
      idx <- idx + 1

      if (verbose && b %% 100 == 0) message(sprintf("MI %d, Bootstrap %d", m, b))
    }

    if (verbose) message(sprintf("Completed MI %d of %d", m, M))
  }

  # NOTE: this uses the (1 - boot_quant) quantile (5% for boot_quant = 0.95), whereas
  # thresh_SI_boot() and thresh_MI() use (1 - boot_quant) / 2 (2.5%). Left unchanged
  # because it alters results; use probs = (1 - boot_quant) / 2 to make them consistent.
  lb_ci <- stats::quantile(thresholds_all, probs = 1 - boot_quant, na.rm = TRUE, names = FALSE)
  flagged_outliers <- which(RD_org > lb_ci)

  result <- list(
    thresholds       = thresholds_all,
    threshold        = lb_ci,
    flagged_outliers = flagged_outliers,
    n_obs            = n_time,
    call             = call
  )

  class(result) <- "MI_boot_result"
  result
}
