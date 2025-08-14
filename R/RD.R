#' Comprehensive Outlier Detection Using Robust Distance Thresholding
#'
#' Performs univariate outlier detection + imputation, robust distance, and multiple thresholding methods.
#'
#' @inheritParams x
#' @inheritParams w
#' @inheritParams method
#' @inheritParams mode
#' @inheritParams cov_mcd
#' @inheritParams ind_incld
#' @inheritParams dist
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
#' @return Depends on method:
#' \describe{
#'   \item{Single method}{Returns the result from the specific threshold method.}
#'   \item{RD_obj}{The robust distance object from compute_RD().}
#'   \item{call}{The matched function call.}
#' }
#'
#' @export
RD <- function(x, w = NULL, method = c("SI_boot", "MI", "MI_boot", "SI","F", "SHASH"),
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
  method <- match.arg(method)

  if (verbose) message("Computing robust distances and covariance.")
  RD_obj <- compute_RD(x = x, mode = mode, cov_mcd = cov_mcd, ind_incld = ind_incld, dist = dist)

  thresholds <- threshold_RD(x = x, w = w, method = method, RD_obj = RD_obj,
                             impute_method = impute_method, cutoff = cutoff, trans = trans,
                             M = M, k = k, alpha = alpha, quantile = quantile,
                             verbose = verbose, boot_quant = boot_quant, B = B)

  if(method == "SHASH"){
    RD_obj = thresholds$RD_obj_shash
  }


  result <- list(
    thresholds = thresholds$thresholds,
    RD_obj = thresholds$RD_obj,
    call = call
  )

  class(result) <- "RD"
  result
}




