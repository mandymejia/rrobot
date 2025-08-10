#' Comprehensive Outlier Detection Using Robust Distance Thresholding
#'
#' Performs univariate outlier detection + imputation, robust distance, and multiple thresholding methods.
#'
#' @param x A numeric matrix or data frame of dimensions T × Q.
#' @param w A numeric matrix (n_time × L) of low-kurtosis predictors (optional).
#' @param threshold_method Character string; one of "all","SI","SI_boot","MI","MI_boot","F".
#' @param mode Character string; either "auto" or "manual" for robust covariance computation.
#' @param cov_mcd Optional covariance matrix (Q × Q); required in "manual" mode.
#' @param ind_incld Optional vector of row indices; required in "manual" mode.
#' @param dist Logical; if TRUE, compute squared robust distances.
#' @param impute_method Character string; imputation method for univariate outliers.
#' @param cutoff Numeric; threshold multiplier for univariate outlier detection.
#' @param trans Character string; transformation method, one of "SHASH" or "robZ".
#' @param M Integer; number of multiple imputation datasets.
#' @param k Integer; number of perturbation cycles per imputation.
#' @param alpha Numeric; significance level for quantile thresholding.
#' @param quantile Numeric; quantile level for F-distribution threshold.
#' @param verbose Logical; if TRUE, print progress messages.
#' @param boot_quant Numeric; confidence level for bootstrap confidence intervals.
#' @param B Integer; number of bootstrap samples.
#'
#' @return Depends on threshold_method:
#' \describe{
#'   \item{Single method}{Returns the result from the specific threshold method.}
#'   \item{"all"}{Returns a list with results from all threshold methods.}
#' }
#'
#' @export
RD <- function(x, w = NULL, threshold_method = c("SI_boot", "MI", "MI_boot", "SI","F"),
                         # RD parameters
                         mode = "auto", cov_mcd = NULL, ind_incld = NULL, dist = TRUE,
                         # impute_univOut paramters
                         impute_method = "mean",
                         # univOut parameters
                         cutoff = 4,
                         trans = "SHASH",
                         # MImpute parameters
                         M = 50,
                         k = 100,
                         # Threshold parameters
                         alpha = 0.01,
                         quantile = 0.01,
                         verbose = FALSE,
                         boot_quant = 0.95,
                         B = 1000) {
  call <- match.call()
  threshold_method <- match.arg(threshold_method)


  RD_obj <- compute_RD(x = x, mode = mode, cov_mcd = cov_mcd, ind_incld = ind_incld, dist = dist)

  thresholds <- threshold_RD(x = x, w = w, threshold_method = threshold_method, RD_obj = RD_obj,
                             impute_method = impute_method, cutoff = cutoff, trans = trans,
                             M = M, k = k, alpha = alpha, quantile = quantile,
                             verbose = verbose, boot_quant = boot_quant, B = B)

  result <- list(
    RD_obj = RD_obj,
    thresholds = thresholds,
    call = call
  )

  class(result) <- "RD"
  result
  }
