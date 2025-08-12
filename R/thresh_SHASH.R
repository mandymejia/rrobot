#' SHASH-based Outlier Detection with F-Distribution Thresholding
#'
#' Performs SHASH-based univariate outlier detection and applies F-distribution thresholding
#' using the SHASH-determined clean subset size.
#'
#' @param x A numeric matrix or data frame of dimensions T × Q.
#' @param cutoff Numeric; threshold multiplier for SHASH-based univariate outlier detection (default = 4).
#' @param quantile Numeric; quantile level for F-distribution threshold (default = 0.01).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{c}{Consistency correction factor for robust distances.}
#'   \item{m}{Estimated degrees of freedom parameter used in the F-distribution.}
#'   \item{df}{A numeric vector of degrees of freedom: \code{c(df1, df2)}.}
#'   \item{scale}{Scale factor for the threshold.}
#'   \item{threshold}{Threshold for squared robust distances.}
#'   \item{call}{The matched function call.}
#' }
#'
#' @details
#' This function implements a simple two-step process:
#' \enumerate{
#'   \item Uses SHASH transformation with isolation forest to detect univariate outliers
#'   \item Applies F-distribution thresholding using SHASH-determined clean subset size
#' }
#'
#' The key difference from regular F-distribution thresholding is that h (the robust subset size)
#' is determined by SHASH outlier detection rather than MCD automatic selection.
#'
#' @seealso \code{\link{univOut}}, \code{\link{thresh_F}}
#' @export
thresh_SASH <- function(x, cutoff = 4, quantile = 0.01) {
  call <- match.call()


  univOut_result <- univOut(x = x, cutoff = cutoff, method = "SHASH")


  imp_result <- impute_univOut(x = x, outlier_mask = univOut_result$outliers, method = "interp")

  RD_obj <- compute_RD(x = imp_result$imp_data, mode = "auto", dist = TRUE)
  F_result <- thresh_F(Q = ncol(x), n = nrow(x), h = RD_obj$h, quantile = quantile)

  result <- list(
    c = F_result$c,
    m = F_result$m,
    df = F_result$df,
    scale = F_result$scale,
    threshold = F_result$threshold,
    call = call
  )

  class(result) <- "SASH_result"
  result
}
