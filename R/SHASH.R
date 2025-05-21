#' Robust empirical rule
#'
#' Robust empirical rule outlier detection applicable to approximately Normal data
#'
#' @param x The data
#' @param thr MAD threshold
#' @param symmetric Single scale for the entire data? Default: \code{TRUE}.
#' @param use_huber Use the Huber estimates for center and scale? Default: 
#'  \code{FALSE}.
#' @param upper_only Only consider upper threshold? Default: \code{FALSE}.
#'
#' @return Logical vector indicating whether each element in \code{x} is an
#'  outlier (\code{TRUE} if an outlier).
#' @importFrom MASS huber
#' @importFrom stats median mad
#' @keywords internal
emprule_rob <- function(x, thr = 4, symmetric = TRUE, use_huber = FALSE, upper_only = FALSE) {
  
  # Validate inputs
  if (!is.numeric(x)) stop("Input data 'x' must be numeric.")
  if (!is.numeric(thr) || length(thr) != 1 || thr <= 0) {
    stop("Threshold 'thr' must be a positive numeric value.")
  }
  
  # Error handling for upper_only mode
  if (use_huber && upper_only) {
    stop("Cannot use `upper_only = TRUE` when `use_huber = TRUE`. Set `use_huber = FALSE`.")
  }
  
  # Error handling for use_huber mode
  if (use_huber && !symmetric) {
    stop("Cannot use `use_huber = TRUE` when `symmetric = FALSE`. Set `symmetric = TRUE`.")
  }
  
  # Calculate the center and scale
  x_med <- stats::median(x, na.rm = TRUE)
  
  if (use_huber) {
    # Use Huber's estimate for location and scale
    huber_fit <- MASS::huber(x)
    center <- huber_fit$mu
    scale <- huber_fit$s
    
  } else if (symmetric) {
    # Use Median Absolute Deviation (MAD) scaled to standard deviation
    MAD <- stats::mad(x, na.rm = TRUE)
    center <- x_med
    scale <- MAD
    
  } else {
    # Asymmetric: Calculate separate left and right scales
    xl <- x[x < x_med]  # Points to the left of the median
    xr <- x[x > x_med]  # Points to the right of the median
    
    left_mad <- 1.4826 * stats::median(abs(xl - x_med), na.rm = TRUE)
    right_mad <- 1.4826 * stats::median(abs(xr - x_med), na.rm = TRUE)
    center <- x_med
  }

  # Calculate thresholds
  if (symmetric) {
    lim_left <- center - thr * scale
    lim_right <- center + thr * scale
  } else {
    lim_left <- center - thr * left_mad
    lim_right <- center + thr * right_mad
  }

  # Identify outliers
  if (upper_only) {
    out <- x > lim_right  # Only filter upper threshold values
  } else {
    out <- (x < lim_left) | (x > lim_right)  # Default: Detect both upper & lower outliers
  }
  
  return(out)
}