#' Outlier Detection via Multiple Imputation Voting (MI)
#'
#' Applies robust distance (RD) computation to multiply imputed datasets, derives thresholds,
#' and flags outliers via majority voting.
#'
#' @param RD_org_obj A list output from `comp_RD()` applied to the original data.
#' @param imp_datasets A list of M numeric matrices (each n_time × Q), multiply imputed versions of original data.
#' @param alpha Significance level used to compute RD threshold (default = 0.01 for 99th percentile).
#'
#' @return A list with:
#' \describe{
#'   \item{thresholds}{A numeric vector of length M; (1 - alpha) quantiles of RD per imputed dataset.}
#'   \item{voted_outliers}{A logical vector of length n_time; TRUE if flagged as outlier in > M/2 imputations.}
#' }
#'
#' @importFrom stats quantile
#' @export
MI <- function(RD_org_obj, imp_datasets, alpha = 0.01) {
  M <- length(imp_datasets)
  n_time <- length(RD_org_obj$RD)
  thresholds <- numeric(M)
  outlier_flags <- matrix(FALSE, nrow = n_time, ncol = M)
  
  cutoff_q <- 1 - alpha
  RD_orig <- RD_org_obj$RD
  cov_mcd <- RD_org_obj$S_star
  ind_incld <- RD_org_obj$ind_incld
  
  for (m in seq_len(M)) {
    imp_data <- imp_datasets[[m]]
    mu_MCD <- colMeans(imp_data[ind_incld, , drop = FALSE])
    xbar_mat <- matrix(mu_MCD, nrow = n_time, ncol = ncol(imp_data), byrow = TRUE)
    invcov_sqrt <- expm::sqrtm(solve(cov_mcd))
    temp <- (imp_data - xbar_mat) %*% invcov_sqrt
    RD_imp <- rowSums(temp^2)
    
    thresholds[m] <- quantile(RD_imp, cutoff_q, na.rm = TRUE)
    outlier_flags[, m] <- RD_orig > thresholds[m]
  }
  
  voted_outliers <- rowSums(outlier_flags) > (M / 2)
  
  return(list(
    thresholds = thresholds,
    voted_outliers = voted_outliers
  ))
}
