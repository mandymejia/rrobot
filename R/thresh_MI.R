#' Outlier Detection via Multiple Imputation Voting (MI)
#'
#' Applies robust distance (RD) computation to multiply imputed datasets, derives thresholds,
#' and flags outliers via majority voting. Also computes the lower bound of the 95% confidence interval
#' of the (1 - alpha) quantiles across imputations.
#'
#' @inheritParams RD_org_obj
#' @inheritParams imp_datasets
#' @inheritParams alpha
#' @inheritParams verbose
#'
#' @return A list with:
#' \describe{
#'   \item{thresholds}{Numeric vector of length M; (1 - alpha) quantiles of RD per imputed dataset.}
#'   \item{LB95_CI}{Lower bound of the 95% confidence interval of thresholds (2.5th percentile).}
#'   \item{call}{The matched function call.}
#' }
#'
#' @importFrom stats quantile
#' @importFrom expm sqrtm
#' @export
thresh_MI <- function(RD_org_obj, imp_datasets, alpha = 0.01, verbose = FALSE) {
  stopifnot(is.list(RD_org_obj), !is.null(RD_org_obj$RD), !is.null(RD_org_obj$S_star), !is.null(RD_org_obj$ind_incld))
  stopifnot(is.list(imp_datasets), length(imp_datasets) > 1)
  if (verbose) message("Running MI method: thresholds across ", length(imp_datasets), " imputations...")

  call <- match.call()

  M <- length(imp_datasets)
  n_time <- length(RD_org_obj$RD)
  p <- ncol(imp_datasets[[1]])

  thresholds <- numeric(M)
  outlier_flags <- matrix(FALSE, nrow = n_time, ncol = M)

  RD_orig <- RD_org_obj$RD
  cov_mcd <- RD_org_obj$S_star
  ind_incld <- RD_org_obj$ind_incld
  cutoff_q <- 1 - alpha

  invcov_sqrt <- RD_org_obj$invcov_sqrt

  for (m in seq_len(M)) {
    imp_data <- imp_datasets[[m]]
    stopifnot(nrow(imp_data) == n_time, ncol(imp_data) == p)

    # mu_MCD <- colMeans(imp_data[ind_incld, , drop = FALSE])
    mu_MCD <- RD_org_obj$xbar_star # Don't re-compute mean
    xbar_mat <- matrix(mu_MCD, nrow = n_time, ncol = p, byrow = TRUE)

    RD_imp <- rowSums(((imp_data - xbar_mat) %*% invcov_sqrt)^2)
    thresholds[m] <- quantile(RD_imp, cutoff_q, na.rm = TRUE)

    outlier_flags[, m] <- RD_orig > thresholds[m]
  }

  LB95_CI <- quantile(thresholds, 0.025, na.rm = TRUE)

  result <- list(
    thresholds = thresholds,
    LB95_CI = LB95_CI,
    call = call
  )

  class(result) <- "MI_result"
  result
}
