#' Outlier Detection via Multiple Imputation Voting (MI)
#'
#' Applies robust distance (RD) computation to multiply imputed datasets, derives thresholds,
#' and flags outliers via majority voting. Also computes the lower bound of the 95% confidence interval
#' of the (1 - alpha) quantiles across imputations.
#'
#' @param RD_org_obj Output list from `comp_RD()` on the original data. Must contain $RD, $S_star, and $ind_incld.
#' @param imp_datasets A list of M numeric matrices (T × Q), multiply imputed versions of original data.
#' @param alpha Significance level used to compute RD threshold (default = 0.01 for 99th percentile).
#'
#' @return A list with:
#' \describe{
#'   \item{thresholds}{Numeric vector of length M; (1 - alpha) quantiles of RD per imputed dataset.}
#'   \item{voted_outliers}{Logical vector (length T); TRUE if RD > threshold in > M/2 imputations.}
#'   \item{LB95_CI}{Lower bound of the 95% confidence interval of thresholds (2.5th percentile).}
#' }
#'
#' @importFrom stats quantile
#' @importFrom expm sqrtm
#' @export
MI <- function(RD_org_obj, imp_datasets, alpha = 0.01) {
  stopifnot(is.list(RD_org_obj), !is.null(RD_org_obj$RD), !is.null(RD_org_obj$S_star), !is.null(RD_org_obj$ind_incld))
  stopifnot(is.list(imp_datasets), length(imp_datasets) > 1)

  M <- length(imp_datasets)
  n_time <- length(RD_org_obj$RD)
  Q <- ncol(imp_datasets[[1]])

  thresholds <- numeric(M)
  outlier_flags <- matrix(FALSE, nrow = n_time, ncol = M)

  RD_orig <- RD_org_obj$RD
  cov_mcd <- RD_org_obj$S_star
  ind_incld <- RD_org_obj$ind_incld
  cutoff_q <- 1 - alpha

  invcov_sqrt <- RD_org_obj$invcov_sqrt

  for (m in seq_len(M)) {
    imp_data <- imp_datasets[[m]]
    stopifnot(nrow(imp_data) == n_time, ncol(imp_data) == Q)

    mu_MCD <- colMeans(imp_data[ind_incld, , drop = FALSE])
    xbar_mat <- matrix(mu_MCD, nrow = n_time, ncol = Q, byrow = TRUE)

    RD_imp <- rowSums(((imp_data - xbar_mat) %*% invcov_sqrt)^2)
    thresholds[m] <- quantile(RD_imp, cutoff_q, na.rm = TRUE)

    outlier_flags[, m] <- RD_orig > thresholds[m]
  }

  voted_outliers <- rowSums(outlier_flags) > (M / 2)
  LB95_CI <- quantile(thresholds, 0.025, na.rm = TRUE)

  result <- list(
    thresholds = thresholds,
    voted_outliers = voted_outliers,
    LB95_CI = LB95_CI
  )

  class(result) <- "MI_result"
  return(result)
}
