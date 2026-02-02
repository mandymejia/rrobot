#' Temporally impute univariate outliers from external detection
#'
#' Takes a high-kurtosis data matrix and a precomputed outlier mask,
#' replaces the outliers with NA, and applies temporal interpolation using `imputeTS::na_interpolation`.
#'
#' @inheritParams x
#' @param outlier_mask A logical matrix (same dimensions as x) with TRUE at outlier positions.
#' @param method One of \code{"mean"} or \code{"interp"}; \code{"interp"} uses \code{imputeTS::na_interpolation},
#'   \code{"mean"} fills NAs with column means.
#'
#' @return A list with elements:
#' \describe{
#'   \item{x_df}{Original matrix with outliers replaced as NA (as tibble).}
#'   \item{NA_data}{Matrix version of x with NAs at outlier positions.}
#'   \item{imp_data}{Imputed matrix after temporal interpolation.}
#'   \item{NA_locs}{Row-column indices of outliers (now NA).}
#'   \item{call}{The matched function call.}
#' }
#'
#' @importFrom tidyr as_tibble
#' @importFrom imputeTS na_interpolation
#' @keywords internal
impute_univOut <- function(x, outlier_mask, method = c("mean", "interp")) {
  stopifnot(is.matrix(x), is.logical(outlier_mask), all(dim(x) == dim(outlier_mask)))

  call <- match.call()
  method <- match.arg(method)

  n_cols <- ncol(x)

  # Step 1: Replace outliers with NA
  NA_x <- x
  NA_x[outlier_mask] <- NA
  NA_x_locs <- which(is.na(NA_x), arr.ind = TRUE)

  # Step 2: Apply temporal interpolation
  if (method == "interp") {
    # Temporal interpolation
    imp_x <- apply(NA_x, 2, imputeTS::na_interpolation)
    imp_x <- as.matrix(imp_x)
  } else if (method == "mean") {
    # compute mean of each column (colMeans with na.rm=TRUE)
    col_means <- colMeans(NA_x, na.rm = TRUE)
    imp_x <- NA_x

    # For each column replace all NAs with that value ^
    for (j in seq_len(n_cols)) {
      na_indices <- is.na(NA_x[, j])
      imp_x[na_indices, j] <- col_means[j]
    }
  }

  # Step 3: Ensure correct dimensions
  if (is.null(dim(imp_x)) || ncol(imp_x) == 1) {
    imp_x <- matrix(imp_x, ncol = n_cols)
  }

  if (!is.numeric(imp_x)) stop("Imputed result contains non-numeric values.")

  # Step 4: Return structured output
  list(
    x_df = suppressMessages(as_tibble(NA_x, .name_repair = "unique")),
    NA_data  = NA_x,
    imp_data = imp_x,
    NA_locs  = NA_x_locs,
    call = call
  )
}
