#' Comprehensive Outlier Detection Using Robust Distance Thresholding
#'
#' Performs univariate outlier detection + imputation, robust distance, and multiple thresholding methods.
#'
#' @param x A numeric matrix or data frame of dimensions T × Q.
#' @param w A numeric matrix (n_time × L) of low-kurtosis predictors (optional).
#' @param threshold_method Character string; one of "all","SI","SI_boot","MI","MI_boot","F", "SASH".
#' @param RD_obj Pre-computed RD object from compute_RD().
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
#' @return A list with:
#' \describe{
#'   \item{thresholds}{Result from the specific threshold method, or list of all methods if "all".}
#'   \item{RD_obj}{The robust distance object from compute_RD().}
#'   \item{call}{The matched function call.}
#' }
#'
#' @export
threshold_RD <- function(x, w = NULL, threshold_method = c("SI_boot", "MI", "MI_boot", "SI","F", "SHASH", "all"), RD_obj = NULL,
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

  # Data pre-processing
  stopifnot("RD_obj is required from compute_RD()" = !is.null(RD_obj))

  out_result <- univOut(x = x, cutoff = cutoff, method = trans)
  imp_result <- impute_univOut(x = x, outlier_mask = out_result$outliers, method = impute_method)


  if (threshold_method %in% c("all", "MI", "MI_boot")) {
    multiple_imp <- MImpute(x = imp_result$imp_data, w = w, outlier_matrix = out_result$outliers, M = M, k = k)
  }

  thresholds <- switch(threshold_method,
                   "SI" = thresh_SI(RD_org_obj = RD_obj, imp_data = imp_result$imp_data, alpha = alpha),

                   "SI_boot" = thresh_SI_boot(RD_org_obj = RD_obj, imp_data = imp_result$imp_data,
                                              B = B, alpha = alpha, boot_quant = boot_quant, verbose = verbose),

                   "MI" = thresh_MI(RD_org_obj = RD_obj, imp_datasets = multiple_imp$imp_datasets, alpha = alpha),

                   "MI_boot" = thresh_MI_boot(RD_org_obj = RD_obj, imp_datasets = multiple_imp$imp_datasets,
                                              B = B, alpha = alpha, boot_quant = boot_quant, verbose = verbose),

                   "F" = thresh_F(Q = ncol(x), n = nrow(x), h = RD_obj$h, quantile = quantile),

                   "SASH" = thresh_SASH(x = x, cutoff = cutoff, quantile = quantile),

                   "all" = list(
                     SI = thresh_SI(RD_org_obj = RD_obj, imp_data = imp_result$imp_data, alpha = alpha),
                     SI_boot = thresh_SI_boot(RD_org_obj = RD_obj, imp_data = imp_result$imp_data,
                                              B = B, alpha = alpha, boot_quant = boot_quant, verbose = verbose),
                     MI = thresh_MI(RD_org_obj = RD_obj, imp_datasets = multiple_imp$imp_datasets, alpha = alpha),
                     MI_boot = thresh_MI_boot(RD_org_obj = RD_obj, imp_datasets = multiple_imp$imp_datasets,
                                              B = B, alpha = alpha, boot_quant = boot_quant, verbose = verbose),
                     F = thresh_F(Q = ncol(x), n = nrow(x), h = RD_obj$h, quantile = quantile),
                     SASH = thresh_SASH(x = x, cutoff = cutoff, quantile = quantile)
                   )
  )


  result <- list(
    thresholds = thresholds,
    RD_obj = RD_obj,
    call = call
  )

  class(result) <- "RD"
  result
}
