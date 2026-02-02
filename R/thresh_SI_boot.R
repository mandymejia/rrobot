#' Compute SI Boot Thresholds for Outlier Detection
#'
#' Computes a robust distance (RD)–based threshold using single-imputed data
#' followed by bootstrap resampling over clean (included) indices. Returns the confidence interval
#' bounds of the bootstrapped 99th percentiles.
#'
#' @inheritParams RD_org_obj
#' @inheritParams imp_data
#' @inheritParams B
#' @inheritParams alpha
#' @inheritParams boot_quant
#' @inheritParams verbose
#'
#' @return A list with:
#' \describe{
#'   \item{thresholds}{Vector of 99th quantiles of RD for each bootstrap sample.}
#'   \item{threshold}{Threshold based on lower bound of the confidence interval.}
#'   \item{flagged_outliers}{Integer vector of row indices from original data matrix that exceed the threshold.}
#'   \item{UB_CI}{Upper bound of the confidence interval for the 99th quantiles.}
#'   \item{call}{The matched function call.}
#' }
#' @keywords internal
thresh_SI_boot <- function(RD_org_obj, imp_data,
                           B = 1000, alpha = 0.01, boot_quant = 0.95,
                           verbose = FALSE) {

  call <- match.call()

  ind_incld <- RD_org_obj$ind_incld
  ind_excld <- RD_org_obj$ind_excld
  center0   <- RD_org_obj$xbar_star
  S0        <- RD_org_obj$S_star

  quant99  <- numeric(B)
  cutoff_q <- 1 - alpha

  if (verbose) message("Running SI_boot method: ", B, " bootstrap samples...")
  for (b in seq_len(B)) {
    # Sample included and excluded with replacement, then combine
    boot_idx_incld <- sample(ind_incld, size = length(ind_incld), replace = TRUE)
    boot_idx_excld <- sample(ind_excld, size = length(ind_excld), replace = TRUE)
    boot_idx <- sort(c(boot_idx_incld, boot_idx_excld))

    imp_boot <- imp_data[boot_idx, , drop = FALSE]

    # Use original robust center/covariance from RD_org_obj (no recomputation)
    RD_boot <- stats::mahalanobis(imp_boot, center = center0, cov = S0)

    quant99[b] <- unname(stats::quantile(RD_boot, cutoff_q, na.rm = TRUE))

    if (verbose && b %% 100 == 0) message("Bootstrap ", b, "/", B, " complete.")
  }

  lower_p <- (1 - boot_quant) / 2
  upper_p <- 1 - lower_p
  LB_CI <- stats::quantile(quant99, probs = lower_p, na.rm = TRUE)
  UB_CI <- stats::quantile(quant99, probs = upper_p, na.rm = TRUE)

  flagged_outliers <- which(RD_org_obj$RD > LB_CI)

  result <- list(
    thresholds   = quant99,
    threshold = LB_CI,
    flagged_outliers = flagged_outliers,
    UB_CI     = UB_CI,
    call      = call
  )
  class(result) <- "SI_boot_result"
  result
}
