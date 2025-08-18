#' Temporal univariate outlier detection using SHASH, robust Yeo-Johnson, or robust MAD.
#'
#' Detects univariate outliers across time for each variable (column) using one of three methods:
#' \itemize{
#'   \item \code{"SHASH"}: Applies SHASH transformation and detects outliers using isolation forest.
#'   \item \code{"robZ"}: Applies robust z-score outlier detection using median and MAD.
#' }
#'
#' @inheritParams x
#' @inheritParams cutoff
#' @inheritParams method_univOut
#'
#' @return A list with elements:
#' \describe{
#'   \item{outliers}{Logical matrix of the same dimension as \code{x}, indicating detected outlier locations (TRUE = outlier).}
#'   \item{method}{A character string indicating the transformation method used.}
#'   \item{call}{The matched function call.}
#' }
#'
#' @export
#' @importFrom stats median
#' @importFrom cellWise transfo
univOut <- function(x, cutoff = 4, method = c("SHASH", "robZ")) {
  call <- match.call()
  method <- match.arg(method)

  n_obv <- nrow(x)
  n_var <- ncol(x)
  outlier_matrix <- matrix(FALSE, n_obv, n_var)
  shash_matrix <- matrix(NA, n_obv, n_var)
  for (ii in seq_len(n_var)) {
    temp <- x[, ii]

    if(method == "SHASH") {
      result <- SHASH_out(x = temp, thr = cutoff, use_isotree = TRUE
                          , thr0 = 3, thr_isotree = 0.6)
      outlier_matrix[result$out_idx, ii] <- TRUE
      shash_matrix[, ii] <- result$x_norm

    }  else if (method == "robZ") {
      temp_med <- median(temp, na.rm = TRUE)
      mad_val <- 1.4826 * median(abs(temp - temp_med), na.rm = TRUE)
      ind_out <- which(temp < (temp_med - cutoff * mad_val) |
                         temp > (temp_med + cutoff * mad_val))
      outlier_matrix[ind_out, ii] <- TRUE
    }
  }

  list(
    # Add function call (params -- include both of them with same argument name)
    outliers = outlier_matrix,
    x_norm = shash_matrix, # normality transformed data matrix when method =SHASH, if method = robZ this returns NA matrix
    method = method,
    call <- call
  )
}
