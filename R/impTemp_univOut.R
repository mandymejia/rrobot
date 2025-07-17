#' Temporally impute outliers for High- and Low-Kurtosis Components
#'
#' Detects univariate outliers in high-kurtosis (hk) and low-kurtosis (lk)
#' components using the robust Yeo–Johnson transformation and a user-specified
#' cutoff, replaces them with NA, then performs temporal interpolation to
#' impute the missing values. Returns both original and imputed data, plus
#' metadata for downstream modeling.
#'
#' @param hk_data Numeric matrix of hk components (T × K).
#' @param lk_data Numeric matrix of lk components (T × L).
#' @param cutoff Numeric threshold for outlier detection
#'   (passed to \code{fMRIscrub:::RD_univOut()}).
#'
#' @return A list with:
#' \describe{
#'   \item{\code{hk_df}}{Original hk data as tibble, columns prefixed `hk_`.}
#'   \item{\code{NA_data}}{hk matrix with outliers replaced by NA.}
#'   \item{\code{imp_data}}{Imputed hk matrix after interpolation.}
#'   \item{\code{comb_data}}{Combined imputed hk & lk tibble, columns `hk_`/`lk_`.}
#'   \item{\code{NA_locs}}{Matrix of (row, col) indices for hk outliers.}
#' }
#'
#' @importFrom fMRIscrub RD_univOut
#' @importFrom imputeTS na_interpolation
#' @importFrom tibble as_tibble
#' @export
impTemp_univOut <- function(hk_data, lk_data, cutoff) {
  # Identify univariate outliers
  univ_hk_out <- RD_univOut(data = hk_data, cutoff = cutoff, trans = "robust-YJ")
  univ_lk_out <- RD_univOut(data = lk_data, cutoff = cutoff, trans = "robust-YJ")

  # Replace TRUE flags with NA
  univ_hk_out[univ_hk_out == TRUE] <- NA
  univ_lk_out[univ_lk_out == TRUE] <- NA

  # Locations of NA entries
  NA_hk_locs <- which(is.na(univ_hk_out), arr.ind = TRUE)
  NA_lk_locs <- which(is.na(univ_lk_out), arr.ind = TRUE)

  # Inject NAs into original data
  hk_data[NA_hk_locs] <- NA
  lk_data[NA_lk_locs] <- NA

  # Perform temporal interpolation
  imp_hk_data <- apply(hk_data, 2, imputeTS::na_interpolation)
  imp_lk_data <- apply(lk_data, 2, imputeTS::na_interpolation)

  # Convert to tibbles and label columns
  hk_data_df <- as_tibble(hk_data)
  hk_imp_df  <- as_tibble(imp_hk_data)
  lk_imp_df  <- as_tibble(imp_lk_data)

  colnames(hk_data_df) <- paste0("hk_", colnames(hk_data_df))
  colnames(hk_imp_df)  <- paste0("hk_", colnames(hk_imp_df))
  colnames(lk_imp_df)  <- paste0("lk_", colnames(lk_imp_df))

  # Combine for downstream modeling
  comb_lk_hk <- cbind(hk_imp_df, lk_imp_df)

  return(list(
    hk_df     = hk_data_df,
    NA_data   = hk_data,
    imp_data  = imp_hk_data,
    comb_data = comb_lk_hk,
    NA_locs   = NA_hk_locs
  ))
}
