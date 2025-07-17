#' Temporal univariate outlier detection using robust MAD.
#'
#' Detects univariate outliers across time for each variable (column) using either no transformation,
#' a robust Yeo-Johnson transformation, or SHASH (not yet implemented).
#'
#' @param data A T × Q matrix of time series (T time points, Q variables).
#' @param cutoff A numeric value indicating how many MADs away from the median to flag as outliers.
#' @param trans Character string. One of \code{"none"}, \code{"robust-YJ"}, or \code{"SHASH"}.
#'
#' @return A logical matrix of the same dimension as \code{data}, indicating outlier locations.
#' @export
#' @importFrom stats median
#' @importFrom cellWise transfo
RD_univOut <- function(data, cutoff = 4, trans = c("none", "robust-YJ", "SHASH")) {
  trans <- match.arg(trans, c("none", "robust-YJ", "SHASH"))
  t <- dim(data)[1]
  Q <- dim(data)[2]
  univOut <- matrix(FALSE, t, Q)

  if (trans == "none") {
    for (ii in seq(Q)) {
      temp <- data[, ii]
      temp_med <- median(temp, na.rm = TRUE)
      mad <- 1.4826 * median(abs(temp - temp_med), na.rm = TRUE)
      cutoff1 <- temp_med - (cutoff * mad)
      cutoff2 <- temp_med + (cutoff * mad)
      ind_out <- which(temp <= cutoff1 | temp >= cutoff2)
      univOut[ind_out, ii] <- TRUE
    }
  } else if (trans == "robust-YJ") {
    for (ii in seq(Q)) {
      temp <- data[, ii]
      trans_temp <- (cellWise::transfo(temp, type = "YJ", robust = TRUE, standardize = TRUE))$Xt
      medi <- median(trans_temp)
      MAD <- median(abs(trans_temp - medi))
      STD <- 1.4826 * MAD
      cutoff1 <- medi - (cutoff * STD)
      cutoff2 <- medi + (cutoff * STD)
      ind_out <- which(trans_temp <= cutoff1 | trans_temp >= cutoff2)
      univOut[ind_out, ii] <- TRUE
    }
  } else {
    stop("SHASH not implemented yet.")
  }

  univOut
}
