#' Temporal univariate outlier detection using SHASH or robust MAD.
#'
#' Detects univariate outliers across time for each variable (column) using one of two methods:
#' \itemize{
#'   \item \code{"SHASH"}: Applies SHASH transformation and detects outliers using isolation forest.
#'   \item \code{"robZ"}: Applies robust z-score outlier detection using median and MAD.
#' }
#'
#' @inheritParams params_doc
#' @inheritParams method_univOut_doc
#'
#' @return A list with elements:
#' \describe{
#'   \item{outliers}{Logical matrix of the same dimension as \code{x}, indicating detected outlier locations (TRUE = outlier).}
#'   \item{x_norm}{Numeric matrix of normality-transformed data when \code{method = "SHASH"}; a matrix of \code{NA} when \code{method = "robZ"}.}
#'   \item{method}{A character string indicating the transformation method used.}
#'   \item{call}{The matched function call.}
#' }
#'
#' @keywords internal
univOut <- function(x, cutoff = 4, method = c("SHASH", "robZ")) {
  call <- match.call()
  method <- match.arg(method)
  x <- as.matrix(x)

  n_obv <- nrow(x)
  n_var <- ncol(x)
  outlier_matrix <- matrix(FALSE, n_obv, n_var)
  shash_matrix <- matrix(NA_real_, n_obv, n_var)
  for (ii in seq_len(n_var)) {
    temp <- x[, ii]

    if (method == "SHASH") {
      result <- SHASH_out(x = temp, thr = cutoff, method_init = "isoplus",
                          thr0 = 3, thr_iso = 0.6)
      outlier_matrix[result$out_idx, ii] <- TRUE
      shash_matrix[, ii] <- result$x_norm

    } else if (method == "robZ") {
      # Same median +/- cutoff * MAD rule as before, but a zero MAD no longer
      # flags every non-median value.
      ind_out <- which(emprule_rob(temp, thr = cutoff, tail = "both"))
      outlier_matrix[ind_out, ii] <- TRUE
    }
  }

  list(
    outliers = outlier_matrix,
    x_norm   = shash_matrix,
    method   = method,
    call     = call
  )
}
