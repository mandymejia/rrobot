#' Comprehensive Outlier Detection Using Robust Distance Thresholding
#'
#' Performs univariate outlier detection + imputation, robust distance, and multiple thresholding methods.
#'
#' @inheritParams params_doc
#' @inheritParams Threshold_method_doc
#'
#' @return A list of class \code{"RD"} with:
#' \describe{
#'   \item{thresholds}{Result from the specific threshold method, or a named list of all methods if \code{"all"}.}
#'   \item{RD_obj}{The robust distance object from compute_RD().}
#'   \item{RD_obj_shash}{Robust distance object computed on the SHASH-normalised data (only for \code{"SHASH"} and \code{"all"}; otherwise NULL).}
#'   \item{x_norm}{SHASH-normalised data matrix (only for \code{"SHASH"} and \code{"all"}; otherwise NULL).}
#'   \item{call}{The matched function call.}
#' }
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(300), ncol = 3)
#' x[1:3, ] <- x[1:3, ] + 6
#' rd <- compute_RD(x)
#' res <- threshold_RD(x, method = "F", RD_obj = rd, trans = "robZ")
#' res$thresholds$threshold
#'
#' @export
threshold_RD <- function(x, w = NULL, method = c("SI_boot", "MI", "MI_boot", "SI", "F", "SHASH", "all"), RD_obj = NULL,
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

  # Input checks (cheap, so they run before any heavy computation)
  stopifnot("When method = 'SHASH', trans must also be 'SHASH'" = !(method == "SHASH" && trans != "SHASH"))
  stopifnot("RD_obj is required from compute_RD()" = !is.null(RD_obj))
  stopifnot("RD_obj must contain squared robust distances; run compute_RD() with dist = TRUE" = !is.null(RD_obj$RD))
  if (method %in% c("all", "MI", "MI_boot") && !is.matrix(w)) {
    stop("method = '", method, "' requires 'w', a numeric matrix of predictors (n_time x L).")
  }

  if (verbose) message("Detecting univariate outliers using ", trans, " method.")
  out_result <- univOut(x = x, cutoff = cutoff, method = trans) # univariate outlier detection
  if (verbose) message("Imputing detected outliers using ", impute_method, " method.")
  imp_result <- impute_univOut(x = x, outlier_mask = out_result$outliers, method = impute_method) # univariate outlier imputation

  if (method %in% c("all", "MI", "MI_boot")) {
    if (verbose) message("Generating ", M, " multiply imputed datasets.")
    multiple_imp <- MImpute(x = imp_result$imp_data, w = w, outlier_matrix = out_result$outliers, M = M, k = k)
  } else {
    multiple_imp <- NULL
  }

  RD_obj_shash <- x_norm <- NULL

  if (method %in% c("all", "SHASH")) {
    if (verbose) message("Computing SHASH-transformed robust distances.")
    x_norm <- out_result$x_norm
    # truncating extreme values to avoid numerical errors
    x_norm[x_norm > 100] <- 100
    x_norm[x_norm < -100] <- -100
    # standardizing the column
    x_norm_scaled <- scale(x_norm)
    RD_obj_shash <- compute_RD(x = x_norm_scaled, mode = "auto", dist = TRUE)
  }

  n <- nrow(x)
  p <- ncol(x)

  run_method <- function(m) {
    switch(m,
      "SI" = thresh_SI(RD_org_obj = RD_obj, imp_data = imp_result$imp_data,
                       alpha = alpha, verbose = verbose),

      "SI_boot" = thresh_SI_boot(RD_org_obj = RD_obj, imp_data = imp_result$imp_data,
                                 B = B, alpha = alpha, boot_quant = boot_quant, verbose = verbose),

      "MI" = thresh_MI(RD_org_obj = RD_obj, imp_datasets = multiple_imp$imp_datasets,
                       alpha = alpha, boot_quant = boot_quant, verbose = verbose),

      "MI_boot" = thresh_MI_boot(RD_org_obj = RD_obj, imp_datasets = multiple_imp$imp_datasets,
                                 B = B, alpha = alpha, boot_quant = boot_quant, verbose = verbose),

      "F" = thresh_F(p = p, n = n, h = RD_obj$h, quantile = quantile,
                     RD_obj = RD_obj, verbose = verbose),

      "SHASH" = thresh_F(p = p, n = n, h = RD_obj_shash$h, quantile = quantile,
                         RD_obj = RD_obj_shash, SHASH = TRUE, verbose = verbose)
    )
  }

  if (method == "all") {
    all_methods <- c("SI", "SI_boot", "MI", "MI_boot", "F", "SHASH")
    thresholds <- stats::setNames(lapply(all_methods, run_method), all_methods)
  } else {
    thresholds <- run_method(method)
  }

  result <- list(
    thresholds = thresholds,  # Shared with RD()
    RD_obj = RD_obj,
    RD_obj_shash = RD_obj_shash,
    x_norm = x_norm,
    call = call
  )

  # Attributes for S3 plotting.
  attr(result, "x_data") <- x
  attr(result, "univOut_hk") <- out_result
  attr(result, "impute_univOut_hk") <- imp_result
  attr(result, "MImpute_hk") <- multiple_imp

  class(result) <- "RD"
  result
}
