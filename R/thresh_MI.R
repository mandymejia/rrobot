#' Outlier Detection via Multiple Imputation Voting (MI)
#'
#' Applies robust distance (RD) computation to multiply imputed datasets, derives thresholds,
#' and flags outliers via majority voting. Also computes the lower bound of the 95% confidence interval
#' of the (1 - alpha) quantiles across imputations.
#'
#' @inheritParams RD_org_obj
#' @inheritParams imp_datasets
#' @inheritParams alpha
#' @inheritParams boot_quant
#' @inheritParams verbose
#'
#' @return A list with:
#' \describe{
#'   \item{thresholds}{Numeric vector of length M; (1 - alpha) quantiles of RD per imputed dataset.}
#'   \item{threshold}{Lower bound of the confidence interval of thresholds.}
#'   \item{call}{The matched function call.}
#'   \item{flagged_outliers}{Integer vector of row indices from original data matrix that exceed the threshold.}
#' }
#'
#' @importFrom stats quantile
#' @importFrom expm sqrtm
#' @export
thresh_MI <- function(RD_org_obj, imp_datasets, alpha = 0.01, boot_quant = 0.95, verbose = FALSE) {
  stopifnot(is.list(RD_org_obj), !is.null(RD_org_obj$RD), !is.null(RD_org_obj$S_star), !is.null(RD_org_obj$ind_incld))
  stopifnot(is.list(imp_datasets), length(imp_datasets) > 1)
  if (verbose) message("Running MI method: thresholds across ", length(imp_datasets), " imputations...")

  call <- match.call()

  M <- length(imp_datasets)
  n_time <- length(RD_org_obj$RD)
  p <- ncol(imp_datasets[[1]])

  thresholds <- numeric(M)
  outlier_flags <- matrix(FALSE, nrow = n_time, ncol = M)

  RD_orig  <- RD_org_obj$RD                 # squared RD from original
  cov_mcd  <- RD_org_obj$S_star             # robust covariance (from covMcd)
  ind_incld <- RD_org_obj$ind_incld
  cutoff_q <- 1 - alpha

  for (m in seq_len(M)) {
    imp_data <- imp_datasets[[m]]
    stopifnot(nrow(imp_data) == n_time, ncol(imp_data) == p)

    # Use original robust mean (no recomputation)
    mu_MCD <- RD_org_obj$xbar_star

    # Squared robust distances for this imputed dataset
    RD_imp <- stats::mahalanobis(imp_data, center = mu_MCD, cov = cov_mcd)

    thresholds[m] <- unname(stats::quantile(RD_imp, cutoff_q, na.rm = TRUE))
    outlier_flags[, m] <- RD_orig > thresholds[m]
  }

  LB_CI <- stats::quantile(thresholds, probs = (1 - boot_quant)/2, na.rm = TRUE)

  flagged_outliers <- which(RD_orig > LB_CI)

  result <- list(
    thresholds = thresholds,
    threshold  = LB_CI,
    flagged_outliers = flagged_outliers,
    call       = call
  )

  class(result) <- "MI_result"
  result
}
