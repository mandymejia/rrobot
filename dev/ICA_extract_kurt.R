#' Extract high and low kurtosis ICA components using fMRIscrub
#'
#' This function performs ICA projection on a multivariate time series matrix using `pscrub()` from the
#' fMRIscrub package and separates high- and low-kurtosis components based on the
#' `highkurt` indicator.
#'
#' @param time_series A numeric matrix of dimension T × V (time × variables), where each column is a univariate time series.
#'
#' @return A list with elements:
#' \describe{
#'   \item{hk}{High-kurtosis ICA component matrix (T × K)}
#'   \item{lk}{Low-kurtosis ICA component matrix (T × L)}
#'   \item{highkurt}{Logical vector indicating which components are high-kurtosis}
#' }
#'
#' @importFrom fMRIscrub pscrub
#' @export
#'
#' @examples
#' \dontrun{
#' result <- ICA_extract_kurt(time_series = fMRIscrub::Dat1)
#' }

ICA_extract_kurt <- function(time_series) {
  stopifnot(is.matrix(time_series))

  # Run pscrub with ICA
  message("Running pscrub ICA projection...")
  ICA_result <- pscrub(time_series, projection = "ICA")

  M <- ICA_result$ICA$M  # Time series matrix (T × num_components)
  highkurt_idx <- ICA_result$ICA$highkurt  # Logical vector

  # Subset components
  hk <- M[, highkurt_idx, drop = FALSE]
  lk <- M[, !highkurt_idx, drop = FALSE]

  return(list(
    hk = hk,
    lk = lk,
    highkurt = highkurt_idx
  ))
}
