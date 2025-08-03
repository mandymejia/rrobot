#' Multiple Imputation for High-Kurtosis ICA Components
#'
#' Performs multiple imputation using perturbed robust regression models.
#'
#' @param x A numeric matrix (n_time × Q) of high-kurtosis ICA components to be imputed.
#' @param w A numeric matrix (n_time × L) of low-kurtosis ICA components used as predictors.
#' @param outlier_matrix A logical matrix (same dim as x) indicating univariate outliers to be imputed.
#' @param M Number of multiply imputed datasets (default = 5).
#' @param k Number of perturbation cycles per imputation (default = 10).
#'
#' @return A list with:
#' \describe{
#'   \item{\code{imp_datasets}}{List of M imputed versions of x}
#'   \item{\code{outlier_matrix}}{Logical matrix of imputed outlier positions}
#' }
#'
#' @importFrom MASS mvrnorm
#' @importFrom robustbase lmrob
#' @importFrom stats coef vcov
#' @export
MImpute <- function(x, w, outlier_matrix, M = 50, k = 100) {
  stopifnot(is.matrix(x), is.matrix(w), is.logical(outlier_matrix))
  stopifnot(all(dim(x) == dim(outlier_matrix)))

  n_time <- nrow(x)
  Q <- ncol(x)
  imp_datasets <- vector("list", M)

  for (m in seq_len(M)) {
    x_imp <- x

    for (j in seq_len(Q)) {
      outlier_idx <- which(outlier_matrix[, j])
      if (length(outlier_idx) == 0) next

      x_j <- x[, j]
      x_j[outlier_idx] <- NA
      observed <- !is.na(x_j)

      if (sum(observed) < 10) next  # Skip if insufficient data

      # Combine predictors: w + other high-kurtosis variables (excluding j)
      if (Q > 1) {
        other_x <- x[, -j, drop = FALSE]
        colnames(other_x) <- paste0("X", setdiff(seq_len(Q), j))
        predictors_all <- cbind(w, other_x)
      } else {
        predictors_all <- w
      }
      colnames(predictors_all) <- paste0("V", seq_len(ncol(predictors_all)))

      # ICA output is already standardized, so no need to rescale again
      predictors_scaled <- predictors_all

      # Use only observed rows for fitting
      df_fit <- data.frame(y = x_j[observed], predictors_scaled[observed, , drop = FALSE])

      fit <- tryCatch(
        robustbase::lmrob(y ~ ., data = df_fit, max.it = 100),
        error = function(e) NULL
      )
      if (is.null(fit)) next

      beta_hat <- coef(fit)
      beta_cov <- tryCatch(vcov(fit), error = function(e) diag(length(beta_hat)) * 1e-4)

      # Create design matrix for outlier rows
      new_X <- cbind(1, predictors_scaled[outlier_idx, , drop = FALSE])
      preds_mat <- matrix(NA, nrow = length(outlier_idx), ncol = k)

      for (iter in seq_len(k)) {
        beta_perturbed <- tryCatch(
          MASS::mvrnorm(1, mu = beta_hat, Sigma = beta_cov),
          error = function(e) beta_hat
        )
        preds_mat[, iter] <- new_X %*% beta_perturbed
      }

      x_imp[outlier_idx, j] <- rowMeans(preds_mat)
    }

    imp_datasets[[m]] <- x_imp
  }

  list(
    imp_datasets = imp_datasets,
    outlier_matrix = outlier_matrix
  )
}
