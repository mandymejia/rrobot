#' Compute SI Threshold for Outlier Detection
#'
#' Computes a robust distance (RD) threshold based on single imputation (SI),
#' using the robust covariance from the original data (via RD_org_obj) and
#' recomputed mean from the imputed data.
#'
#' @inheritParams RD_org_obj
#' @inheritParams imp_data
#' @inheritParams alpha
#'
#' @return A list with:
#' \describe{
#'   \item{SI_obj}{A list from \code{\link{compute_RD}} containing robust distances.}
#'   \item{SI_threshold}{Numeric threshold based on the (1 - alpha) quantile of RD.}
#'   \item{call}{The matched function call.}
#' }
#' @export
thresh_SI <- function(RD_org_obj, imp_data, alpha = 0.01) {
  call <- match.call()

  cov_mcd <- RD_org_obj$S_star
  ind_incld <- RD_org_obj$ind_incld

  # Compute robust distances using fixed covariance, recomputed mean
  SI_obj <- compute_RD(
    x = imp_data,
    mode = "manual",
    cov_mcd = cov_mcd,
    ind_incld = ind_incld
  )

  # Threshold = (1 - alpha) quantile of RD
  SI_threshold <- quantile(SI_obj$RD, 1 - alpha, na.rm = TRUE)

  result <- list(
    SI_obj = SI_obj,
    SI_threshold = SI_threshold,
    call = call
  )

  class(result) <- "SI_result"
  result
}
