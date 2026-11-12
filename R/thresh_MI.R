#' Outlier Detection via Multiple Imputation (MI)
#'
#' Applies robust distance (RD) computation to multiply imputed datasets and
#' derives a threshold from the lower bound of the confidence interval of the
#' (1 - alpha) quantiles across imputations. Outliers are the observations whose
#' original RD exceeds that threshold.
#'
#' @inheritParams params_doc
#'
#' @return A list of class \code{"MI_result"} with:
#' \describe{
#'   \item{thresholds}{Numeric vector of length M; (1 - alpha) quantiles of RD per imputed dataset.}
#'   \item{threshold}{Lower bound of the confidence interval of thresholds.}
#'   \item{flagged_outliers}{Integer vector of row indices from original data matrix that exceed the threshold.}
#'   \item{call}{The matched function call.}
#' }
#'
#' @importFrom stats quantile mahalanobis
#' @keywords internal
thresh_MI <- function(RD_org_obj, imp_datasets, alpha = 0.01, boot_quant = 0.95, verbose = FALSE) {
  call <- match.call()
  stopifnot(is.list(RD_org_obj), !is.null(RD_org_obj$RD), !is.null(RD_org_obj$S_star), !is.null(RD_org_obj$ind_incld))
  stopifnot(is.list(imp_datasets), length(imp_datasets) > 1)
  if (verbose) message("Running MI method: thresholds across ", length(imp_datasets), " imputations...")

  M <- length(imp_datasets)
  n_time <- length(RD_org_obj$RD)
  p <- ncol(imp_datasets[[1]])

  RD_orig  <- RD_org_obj$RD                 # squared RD from original
  cov_mcd  <- RD_org_obj$S_star             # robust covariance (from covMcd)
  mu_MCD   <- RD_org_obj$xbar_star          # original robust mean (no recomputation)
  cutoff_q <- 1 - alpha

  thresholds <- vapply(imp_datasets, function(imp_data) {
    stopifnot(nrow(imp_data) == n_time, ncol(imp_data) == p)
    # Squared robust distances for this imputed dataset
    RD_imp <- stats::mahalanobis(imp_data, center = mu_MCD, cov = cov_mcd)
    stats::quantile(RD_imp, cutoff_q, na.rm = TRUE, names = FALSE)
  }, numeric(1))

  LB_CI <- stats::quantile(thresholds, probs = (1 - boot_quant) / 2, na.rm = TRUE, names = FALSE)

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
