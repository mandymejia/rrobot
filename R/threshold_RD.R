#' Comprehensive Outlier Detection Using Robust Distance Thresholding
#'
#' Performs univariate outlier detection + imputation, robust distance, and multiple thresholding methods.
#'
#' @inheritParams x
#' @inheritParams w
#' @inheritParams threshold_method
#' @inheritParams RD_obj
#' @inheritParams impute_method
#' @inheritParams cutoff
#' @inheritParams trans
#' @inheritParams M
#' @inheritParams k
#' @inheritParams alpha
#' @inheritParams quantile
#' @inheritParams verbose
#' @inheritParams boot_quant
#' @inheritParams B
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
  stopifnot("When threshold_method = 'SHASH', trans must also be 'SHASH'" = !(threshold_method == "SHASH" && trans != "SHASH"))
  stopifnot("RD_obj is required from compute_RD()" = !is.null(RD_obj))

  out_result <- univOut(x = x, cutoff = cutoff, method = trans) # univariate outlier detection
  imp_result <- impute_univOut(x = x, outlier_mask = out_result$outliers, method = impute_method) # univariate outlier imputation


  if (threshold_method %in% c("all", "MI", "MI_boot")) {
    multiple_imp <- MImpute(x = imp_result$imp_data, w = w, outlier_matrix = out_result$outliers, M = M, k = k)
  }

  RD_obj_shash <- NULL
  if (threshold_method %in% c("all", "SHASH")) {
    x_norm <- out_result$x_norm
    # truncating extreme values to avoid numerical errors
    x_norm[ x_norm > 100] <- 100
    x_norm[ x_norm < -100] <- -100
    RD_obj_shash <- compute_RD(x = x_norm, mode = "auto", dist = TRUE)
  }

  thresholds <- switch(threshold_method,
                   "SI" = thresh_SI(RD_org_obj = RD_obj, imp_data = imp_result$imp_data, alpha = alpha),

                   "SI_boot" = thresh_SI_boot(RD_org_obj = RD_obj, imp_data = imp_result$imp_data,
                                              B = B, alpha = alpha, boot_quant = boot_quant, verbose = verbose),

                   "MI" = thresh_MI(RD_org_obj = RD_obj, imp_datasets = multiple_imp$imp_datasets, alpha = alpha),

                   "MI_boot" = thresh_MI_boot(RD_org_obj = RD_obj, imp_datasets = multiple_imp$imp_datasets,
                                              B = B, alpha = alpha, boot_quant = boot_quant, verbose = verbose),

                   "F" = thresh_F(p = ncol(x), n = nrow(x), h = RD_obj$h, quantile = quantile),

                   "SHASH" = thresh_F(p = ncol(x), n = nrow(x), h = RD_obj_shash$h, quantile = quantile),

                   "all" = list(
                     SI = thresh_SI(RD_org_obj = RD_obj, imp_data = imp_result$imp_data, alpha = alpha),
                     SI_boot = thresh_SI_boot(RD_org_obj = RD_obj, imp_data = imp_result$imp_data,
                                              B = B, alpha = alpha, boot_quant = boot_quant, verbose = verbose),
                     MI = thresh_MI(RD_org_obj = RD_obj, imp_datasets = multiple_imp$imp_datasets, alpha = alpha),
                     MI_boot = thresh_MI_boot(RD_org_obj = RD_obj, imp_datasets = multiple_imp$imp_datasets,
                                              B = B, alpha = alpha, boot_quant = boot_quant, verbose = verbose),
                     F = thresh_F(p = ncol(x), n = nrow(x), h = RD_obj$h, quantile = quantile),
                     SHASH = thresh_F(p = ncol(x), n = nrow(x), h = RD_obj_shash$h, quantile = quantile)
                   )
  )

  result <- list(
    thresholds = thresholds,
    RD_obj = RD_obj,
    RD_obj_shash = RD_obj_shash,
    call = call
  )

  class(result) <- "RD"
  result
}
