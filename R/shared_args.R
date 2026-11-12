#' Shared parameter documentation
#'
#' Internal; used only via \code{@inheritParams params_doc}.
#'
#' @param ... Additional arguments to method-specific functions.
#' @param alpha Significance level used to compute RD threshold (default = 0.01 for 99th percentile).
#' @param B Integer; number of bootstrap samples (per imputed dataset for the MI methods; default = 1000).
#' @param boot_quant Numeric; confidence level for bootstrap confidence intervals (default = 0.95, for 95% CI).
#' @param cov_mcd Optional covariance matrix (p x p); required in "manual" mode.
#' @param cutoff A numeric value indicating how many MADs away from the median to flag as outliers. The default value is set to be 4.
#' @param dist Logical; if TRUE, compute squared robust Mahalanobis distances for all observations.
#' @param imp_data A numeric matrix (T x p) of single-imputed data.
#' @param imp_datasets A list of M numeric matrices (T x p); multiply imputed datasets.
#' @param impute_method Character string; imputation method for univariate outliers.
#' @param ind_incld Optional vector of row indices used to compute the robust mean; required in "manual" mode.
#' @param k Integer; number of chained-equation cycles per imputed dataset (default = 100 in threshold_RD()).
#' @param quantile Numeric in (0,1) specifying the upper quantile for thresholding; the expected False Positive Rate for the chosen threshold.
#' @param M Integer; number of multiply imputed datasets (default = 50 in threshold_RD()).
#' @param mode Character string; either "auto" (default) to compute MCD internally or "manual" to use user-supplied values.
#' @param RD_obj Pre-computed RD_result object from \code{\link{compute_RD}}.
#' @param RD_org_obj Output list from \code{\link{compute_RD}} on the original data. Must contain $RD, $S_star, and $ind_incld.
#' @param thr Threshold multiplier for outlier detection (default = 4).
#' @param trans Character string; transformation method, one of "SHASH" or "robZ".
#' @param verbose Logical; if TRUE, print progress messages.
#' @param w A numeric matrix (n_time x L) of low-kurtosis ICA components used as predictors (required for MI).
#' @param x A numeric matrix or data frame of dimensions T x p (observations x variables).
#' @param thresh_result A threshold result object from any threshold method containing threshold information.
#' @param binwidth Histogram bin width (default = 0.1).
#' @keywords internal
#' @name params_doc
params_doc <- function(..., alpha, B, boot_quant, cov_mcd, cutoff, dist,
                       imp_data, imp_datasets, impute_method, ind_incld, k,
                       quantile, M, mode, RD_obj, RD_org_obj, thr,
                       trans, verbose, w, x, thresh_result, binwidth) NULL


#' Method_univOut parameter documentation
#'
#' @param method Character string. One of \code{"SHASH"} or \code{"robZ"}.
#'
#' @keywords internal
#' @name method_univOut_doc
method_univOut_doc <- function(method) NULL

#' Threshold_method parameter documentation
#'
#' @param method Character string; one of "all","SI","SI_boot","MI","MI_boot","F", "SHASH".
#'
#' @keywords internal
#' @name Threshold_method_doc
Threshold_method_doc <- function(method) NULL
