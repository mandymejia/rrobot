#' Compute SI Boot Thresholds for Outlier Detection
#'
#' Computes a robust distance (RD)-based threshold using single-imputed data
#' followed by bootstrap resampling. Included and excluded indices (from the
#' robust subset) are resampled separately, with replacement. The lower bound of the
#' confidence interval of the bootstrapped (1 - alpha) quantiles is the threshold.
#'
#' @inheritParams params_doc
#'
#' @return A list of class \code{"SI_boot_result"} with:
#' \describe{
#'   \item{thresholds}{Vector of (1 - alpha) quantiles of RD for each bootstrap sample.}
#'   \item{threshold}{Threshold based on the lower bound of the confidence interval.}
#'   \item{flagged_outliers}{Integer vector of row indices from original data matrix that exceed the threshold.}
#'   \item{UB_CI}{Upper bound of the confidence interval for the (1 - alpha) quantiles.}
#'   \item{call}{The matched function call.}
#' }
#' @importFrom stats quantile mahalanobis
#' @keywords internal
thresh_SI_boot <- function(RD_org_obj, imp_data,
                           B = 1000, alpha = 0.01, boot_quant = 0.95,
                           verbose = FALSE) {

  call <- match.call()
  stopifnot(!is.null(RD_org_obj$RD), !is.null(RD_org_obj$S_star))

  ind_incld <- RD_org_obj$ind_incld
  ind_excld <- RD_org_obj$ind_excld
  center0   <- RD_org_obj$xbar_star
  S0        <- RD_org_obj$S_star

  n_incld  <- length(ind_incld)
  n_excld  <- length(ind_excld)
  cutoff_q <- 1 - alpha

  # Distances are computed row by row, so computing them once and resampling
  # the resulting vector is identical to recomputing them for every bootstrap
  # sample (and much faster). Original robust center/covariance, no recomputation.
  RD_imp <- stats::mahalanobis(imp_data, center = center0, cov = S0)

  quant99 <- numeric(B)

  if (verbose) message("Running SI_boot method: ", B, " bootstrap samples...")
  for (b in seq_len(B)) {
    # Sample included and excluded indices with replacement, then combine.
    # (sample.int avoids sample()'s 1:n behaviour when an index set has length 1)
    boot_idx <- c(ind_incld[sample.int(n_incld, replace = TRUE)],
                  ind_excld[sample.int(n_excld, replace = TRUE)])

    quant99[b] <- stats::quantile(RD_imp[boot_idx], cutoff_q, na.rm = TRUE)

    if (verbose && b %% 100 == 0) message("Bootstrap ", b, "/", B, " complete.")
  }

  lower_p <- (1 - boot_quant) / 2
  upper_p <- 1 - lower_p
  LB_CI <- stats::quantile(quant99, probs = lower_p, na.rm = TRUE)
  UB_CI <- stats::quantile(quant99, probs = upper_p, na.rm = TRUE)

  flagged_outliers <- which(RD_org_obj$RD > LB_CI)

  result <- list(
    thresholds       = quant99,
    threshold        = LB_CI,
    flagged_outliers = flagged_outliers,
    UB_CI            = UB_CI,
    call             = call
  )
  class(result) <- "SI_boot_result"
  result
}
