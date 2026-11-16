# -----------------------------------------------------------------------------
# Robust Empirical Rule
# -----------------------------------------------------------------------------

#' Robust Empirical Rule Outlier Detection
#'
#' Detects outliers using the median +/- \code{thr} x MAD rule, where MAD is
#' normalised by 1.4826 to be consistent with the standard deviation under
#' normality.
#'
#' @param x Numeric vector.
#' @param thr Positive numeric scalar. Threshold multiplier for the MAD rule
#'   (default: 4).
#' @param tail Character string: one of \code{"both"} (default),
#'   \code{"upper"}, or \code{"lower"}, indicating which tail(s) to flag.
#'
#' @return A logical vector the same length as \code{x}. \code{TRUE} indicates
#'   an outlier, \code{FALSE} indicates an inlier.
#'
#' @importFrom stats median
#' @keywords internal
emprule_rob <- function(x, thr = 4, tail = c("both", "upper", "lower")) {

  tail <- match.arg(tail)

  x_med  <- median(x, na.rm = TRUE)
  mad_val <- 1.4826 * median(abs(x - x_med), na.rm = TRUE)

  # If MAD is zero (e.g. more than half the data are identical), the rule
  # cannot distinguish outliers - return a safe all-FALSE vector.
  if (mad_val == 0 || is.na(mad_val)) {
    return(rep(FALSE, length(x)))
  }

  z <- (x - x_med) / mad_val

  if (tail == "upper") {
    out <- z > thr
  } else if (tail == "lower") {
    out <- z < -thr
  } else {
    out <- abs(z) > thr
  }

  out
}
