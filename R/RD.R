#' Comprehensive Outlier Detection Using Robust Distance Thresholding
#'
#' Performs univariate outlier detection + imputation, robust distance, and multiple thresholding methods.
#'
#' @inheritParams params_doc
#' @inheritParams Threshold_method_doc
#'
#' @return A list of class \code{"RD"} with:
#' \describe{
#'   \item{thresholds}{The result from the chosen threshold method.}
#'   \item{RD_obj}{The robust distance object from compute_RD().}
#'   \item{RD_obj_shash}{Robust distances of the SHASH-normalised data (only for \code{method = "SHASH"}; otherwise NULL).}
#'   \item{outliers}{Logical vector indicating which observations have RD greater than the threshold.
#'     For \code{method = "SHASH"} the distances are those of the SHASH-normalised data.}
#'   \item{call}{The matched function call.}
#' }
#'
#' @keywords internal
RD <- function(x, w = NULL, method = c("SI_boot", "MI", "MI_boot", "SI", "F", "SHASH"),
               # RD parameters
               mode = "auto", cov_mcd = NULL, ind_incld = NULL, dist = TRUE,
               # impute_univOut parameters
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
  method <- match.arg(method)
  x <- as.matrix(x)

  if (verbose) message("Computing robust distances and covariance.")
  RD_obj <- compute_RD(x = x, mode = mode, cov_mcd = cov_mcd, ind_incld = ind_incld, dist = dist)

  thresholds <- threshold_RD(x = x, w = w, method = method, RD_obj = RD_obj,
                             impute_method = impute_method, cutoff = cutoff, trans = trans,
                             M = M, k = k, alpha = alpha, quantile = quantile,
                             verbose = verbose, boot_quant = boot_quant, B = B)

  # The SHASH threshold lives on the SHASH-normalised distances.
  RD_used <- if (method == "SHASH") thresholds$RD_obj_shash else thresholds$RD_obj
  outliers <- RD_used$RD > thresholds$thresholds$threshold

  result <- list(
    thresholds = thresholds$thresholds, # Shared with threshold_RD()
    RD_obj = thresholds$RD_obj,
    RD_obj_shash = thresholds$RD_obj_shash,
    outliers = outliers,
    call = call
  )

  # Attributes for S3 plotting.
  attr(result, "x_data") <- attr(thresholds, "x_data")
  attr(result, "univOut_hk") <- attr(thresholds, "univOut_hk")
  attr(result, "impute_univOut_hk") <- attr(thresholds, "impute_univOut_hk")
  attr(result, "MImpute_hk") <- attr(thresholds, "MImpute_hk")

  class(result) <- "RD"
  result
}
