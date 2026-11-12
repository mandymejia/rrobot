#' Compute SI Threshold for Outlier Detection
#'
#' Computes a robust distance (RD) threshold based on single imputation (SI),
#' using the robust covariance from the original data (via RD_org_obj) and
#' recomputed mean from the imputed data.
#'
#' @inheritParams params_doc
#'
#' @return A list of class \code{"SI_result"} with:
#' \describe{
#'   \item{SI_obj}{A list from \code{\link{compute_RD}} containing robust distances.}
#'   \item{threshold}{Numeric threshold based on the (1 - alpha) quantile of RD.}
#'   \item{flagged_outliers}{Integer vector of row indices from original data matrix that exceed the threshold.}
#'   \item{call}{The matched function call.}
#' }
#' @importFrom stats quantile
#' @keywords internal
thresh_SI <- function(RD_org_obj, imp_data, alpha = 0.01, verbose = FALSE) {
  call <- match.call()
  if (verbose) message("Running SI method: computing single imputation threshold...")

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
  threshold <- stats::quantile(SI_obj$RD, 1 - alpha, na.rm = TRUE)

  flagged_outliers <- which(RD_org_obj$RD > threshold)

  result <- list(
    SI_obj = SI_obj,
    threshold = threshold,
    flagged_outliers = flagged_outliers,
    call = call
  )

  class(result) <- "SI_result"
  result
}
