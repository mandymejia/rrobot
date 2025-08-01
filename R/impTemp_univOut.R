#' Temporally impute univariate outliers from external detection
#'
#' Takes a high-kurtosis data matrix and a precomputed outlier mask,
#' replaces the outliers with NA, and applies temporal interpolation using `imputeTS::na_interpolation`.
#'
#' @param x A numeric matrix of dimension n_time × n_var (n_time time points, n_var variables).
#' @param outlier_mask A logical matrix (same dimensions as x) with TRUE at outlier positions.
#'
#' @return A list with elements:
#' \describe{
#'   \item{x_df}{Original matrix with outliers replaced as NA (as tibble).}
#'   \item{NA_data}{Matrix version of x with NAs at outlier positions.}
#'   \item{imp_data}{Imputed matrix after temporal interpolation.}
#'   \item{NA_locs}{Row-column indices of outliers (now NA).}
#' }
#'
#' @importFrom tidyr as_tibble
#' @importFrom imputeTS na_interpolation
#' @export
impTemp_univOut <- function(x, outlier_mask) {
  stopifnot(is.matrix(x), is.logical(outlier_mask), all(dim(x) == dim(outlier_mask)))
  
  # Step 1: Replace outliers with NA
  NA_x <- x
  NA_x[outlier_mask] <- NA
  NA_x_locs <- which(is.na(NA_x), arr.ind = TRUE)
  
  # Step 2: Apply temporal interpolation
  imp_x <- apply(NA_x, 2, imputeTS::na_interpolation)
  imp_x <- as.matrix(imp_x)
  
  # Step 3: Ensure correct dimensions
  if (is.null(dim(imp_x)) || ncol(imp_x) == 1) {
    imp_x <- matrix(imp_x, ncol = 1)
  }
  
  if (!is.numeric(imp_x)) stop("Imputed result contains non-numeric values.")
  
  # Step 4: Return structured output
  return(list(
    x_df     = as_tibble(NA_x, .name_repair = "unique"),
    NA_data  = NA_x,
    imp_data = imp_x,
    NA_locs  = NA_x_locs
  ))
}
