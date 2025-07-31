#' Multiple Imputation for High-Kurtosis ICA Components
#'
#' Performs multiple imputation using perturbed robust regression models.
#' For each variable (column) with flagged univariate outliers, performs M imputations
#' using robust linear regression with low-kurtosis components and other high-kurtosis
#' variables as predictors. Adds Gaussian noise to regression coefficients to introduce
#' variation across imputations.
#'
#' @param x A numeric matrix (n_time × Q) of high-kurtosis components to be imputed.
#' @param w A numeric matrix (n_time × L) of low-kurtosis components used as predictors.
#' @param outlier_matrix A logical matrix of the same dimension as \code{x}, indicating
#'   the location of univariate outliers to be imputed.
#' @param M Integer; number of multiple imputations (default = 5).
#' @param k Integer; number of MCMC-like perturbation cycles per imputation (default = 10).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{imp_datasets}}{A list of \code{M} numeric matrices (each of dimension \code{n_time × Q}) representing multiply imputed versions of \code{x}. In each matrix, previously flagged outliers have been replaced by robust regression-based imputations.}
#'   \item{\code{outlier_matrix}}{A logical matrix (same dimension as \code{x}) indicating the positions of univariate outliers that were imputed. This is the same as the input \code{outlier_matrix}, returned for reference.}
#' }
#'
#' @details
#' For each variable \code{x[, j]} with outliers, the function:
#' \enumerate{
#'   \item Replaces outliers with NA in \code{x[, j]}.
#'   \item Uses the remaining variables (\code{w} followed by \code{x[, -j]}) as predictors in a robust regression.
#'   \item Fits a linear model and perturbs the coefficients with Gaussian noise for \code{M × k} cycles.
#'   \item Averages imputed values over \code{k} iterations to stabilize each imputation.
#' }
#'
#' @importFrom MASS mvrnorm
#' @importFrom robustbase lmrob
#' @export
MImpute <- function(x, w, outlier_matrix, M = 5, k = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
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
      
      x_na <- x[, j]
      x_na[outlier_idx] <- NA
      
      predictors <- cbind(w, x[, -j, drop = FALSE])  # Ensure w is prioritized
      predictors <- scale(predictors)
      observed <- !is.na(x_na)
      
      # Fit robust regression model
      fit <- tryCatch(
        robustbase::lmrob(x_na ~ predictors, max.it = 100),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      
      beta_hat <- coef(fit)
      beta_cov <- vcov(fit)
      
      preds_mat <- matrix(NA, nrow = length(outlier_idx), ncol = k)
      
      for (iter in seq_len(k)) {
        beta_perturbed <- tryCatch(
          MASS::mvrnorm(1, mu = beta_hat, Sigma = beta_cov),
          error = function(e) beta_hat
        )
        new_X <- cbind(1, predictors[outlier_idx, , drop = FALSE])
        preds_mat[, iter] <- new_X %*% beta_perturbed
      }
      
      x_imp[outlier_idx, j] <- rowMeans(preds_mat)
    }
    
    imp_datasets[[m]] <- x_imp
  }
  
  return(list(
    imp_datasets = imp_datasets,
    outlier_matrix = outlier_matrix
  ))
}
