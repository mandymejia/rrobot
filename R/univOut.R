#' Temporal univariate outlier detection using SHASH, robust Yeo-Johnson, or robust MAD.
#'
#' Detects univariate outliers across time for each variable (column) using one of three methods:
#' \itemize{
#'   \item \code{"SHASH"}: Applies SHASH transformation and detects outliers using isolation forest.
#'   \item \code{"robust-YJ"}: Applies a robust Yeo-Johnson transformation, then flags based on MAD.
#'   \item \code{"robMAD"}: Applies robust MAD outlier detection without transformation.
#' }
#'
#' @param hk_data A n_time × n_var matrix of high kurtosis (n_time number of time points, n_var number of variables).
#' @param cutoff A numeric value indicating how many MADs away from the median to flag as outliers. The default value is set to be 4.
#' @param trans Character string. One of \code{"SHASH"}, \code{"robust-YJ"}, or \code{"robMAD"}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{outliers}{Logical matrix of the same dimension as \code{hk_data}, indicating detected outlier locations (TRUE = outlier).}
#'   \item{method}{A character string indicating the transformation method used.}
#' }
#'
#' @export
#' @importFrom stats median
#' @importFrom cellWise transfo
univOut <- function(hk_data, cutoff = 4, trans = c("SHASH", "robust-YJ", "robMAD")) {
  trans <- match.arg(trans)
  n_time <- nrow(hk_data)
  n_var <- ncol(hk_data)
  outlier_matrix <- matrix(FALSE, n_time, n_var)
  
  for (ii in seq_len(n_var)) {
    temp <- hk_data[, ii]
    
    if(trans == "SHASH") {
      result <- SHASH_out(x = temp, thr = cutoff, use_isotree = TRUE)
      outlier_matrix[result$out_idx, ii] <- TRUE
      
    } else if (trans == "robust-YJ") {
      trans_temp <- cellWise::transfo(temp, type = "YJ", robust = TRUE, standardize = TRUE)$Y
      temp_med <- median(trans_temp, na.rm = TRUE)
      mad_val <- 1.4826 * median(abs(trans_temp - temp_med), na.rm = TRUE)
      ind_out <- which(trans_temp < (temp_med - cutoff * mad_val) |
                         trans_temp > (temp_med + cutoff * mad_val))
      outlier_matrix[ind_out, ii] <- TRUE
      
    } else if (trans == "robMAD") {
      temp_med <- median(temp, na.rm = TRUE)
      mad_val <- 1.4826 * median(abs(temp - temp_med), na.rm = TRUE)
      ind_out <- which(temp < (temp_med - cutoff * mad_val) |
                         temp > (temp_med + cutoff * mad_val))
      outlier_matrix[ind_out, ii] <- TRUE
    }
  }

  list(
    outliers = outlier_matrix,
    method = trans
  )
}
