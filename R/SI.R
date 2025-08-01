#' Compute SI Threshold for Outlier Detection
#'
#' Computes a robust distance (RD) threshold based on single imputation (SI),
#' using the robust covariance from the original data (via RD_org_obj) and
#' recomputed mean from the imputed data.
#'
#' @param RD_org_obj A list from `comp_RD()` containing original RD, S_star, ind_incld.
#' @param imp_data A numeric matrix (T × Q) of imputed data.
#' @param alpha Significance level for thresholding (e.g., 0.01 for 99th percentile).
#'
#' @return A list with:
#' \describe{
#'   \item{SI_obj}{A list from \code{comp_RD()} containing robust distances.}
#'   \item{SI_threshold}{Numeric threshold based on the (1 - alpha) quantile of RD.}
#' }
#' @export
SI <- function(RD_org_obj, imp_data, alpha = 0.01) {
  cov_mcd <- RD_org_obj$S_star
  ind_incld <- RD_org_obj$ind_incld
  
  # Compute robust distances using fixed covariance, recomputed mean
  SI_obj <- comp_RD(
    data_matrix = imp_data,
    mode = "manual",
    cov_mcd = cov_mcd,
    ind_incld = ind_incld
  )
  
  # Threshold = (1 - alpha) quantile of RD
  SI_threshold <- quantile(SI_obj$RD, 1 - alpha, na.rm = TRUE)
  
  return(list(
    SI_obj = SI_obj,
    SI_threshold = SI_threshold
  ))
}
