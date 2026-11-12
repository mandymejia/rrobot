# File: thresh_F.R

#' Fit F-distribution Parameters for MCD-based Robust Distances
#'
#' Computes the scaling constant and degrees of freedom for the F-distribution
#' approximation of squared robust Mahalanobis distances based on the
#' Minimum Covariance Determinant (MCD) estimator, following the method
#' of Hardin & Rocke (2005).
#'
#' This function is useful for deriving robust outlier detection thresholds
#' in high-dimensional multivariate data contaminated by outliers.
#'
#' @param p Integer. The number of variables (dimension of the data).
#' @param n Integer. The total sample size.
#' @param h Integer. The number of observations retained in the MCD subset.
#' @param SHASH Boolean. If running SHASH variant.
#' @inheritParams params_doc
#'
#' @return A list of class \code{"F_result"} (or \code{"SHASH_result"} when
#' \code{SHASH = TRUE}) with the following elements:
#' \describe{
#'   \item{c}{Consistency correction factor for robust distances.}
#'   \item{m}{Estimated degrees of freedom parameter used in the F-distribution.}
#'   \item{df}{A numeric vector of degrees of freedom: \code{c(df1, df2)}.}
#'   \item{scale}{Scale factor for the threshold.}
#'   \item{threshold}{Threshold for squared robust distances.}
#'   \item{flagged_outliers}{Integer vector of row indices from original data matrix that exceed the threshold.}
#'   \item{call}{The matched function call.}
#' }
#'
#' @importFrom stats pchisq qchisq qf
#' @keywords internal
thresh_F <- function(p, n, h, quantile, RD_obj, SHASH = FALSE, verbose = FALSE) {
  call <- match.call()
  if (verbose) message("Running ", if (SHASH) "SHASH" else "F", " method: F-distribution threshold...")

  # Step 1: consistency correction factor
  c_corr <- stats::pchisq(q = stats::qchisq(p = h / n, df = p), df = p + 2) / (h / n)

  # Step 2: asymptotic approximation for degrees of freedom
  alpha_h <- (n - h) / n
  q_alpha <- stats::qchisq(p = 1 - alpha_h, df = p)

  c_alpha <- (1 - alpha_h) / stats::pchisq(q = q_alpha, df = p + 2)
  c2 <- -0.5 * stats::pchisq(q = q_alpha, df = p + 2)
  c3 <- -0.5 * stats::pchisq(q = q_alpha, df = p + 4)
  c4 <- 3 * c3

  b1 <- c_alpha * (c3 - c4) / (1 - alpha_h)
  b2 <- 0.5 + c_alpha * (c3 - (q_alpha / p) * (c2 + (1 - alpha_h) / 2)) / (1 - alpha_h)

  v1 <- (1 - alpha_h) * b1^2 * (alpha_h * (c_alpha * q_alpha / p - 1)^2 - 1) -
    2 * c3 * c_alpha^2 * (3 * (b1 - p * b2)^2 + (p + 2) * b2 * (2 * b1 - p * b2))

  v2 <- n * (b1 * (b1 - p * b2) * (1 - alpha_h))^2 * c_alpha^2
  v <- v1 / v2

  m <- 2 / (c_alpha^2 * v)

  # Step 3: finite-sample correction
  m <- m * exp(0.725 - 0.00663 * p - 0.078 * log(n))

  df_F <- c(p, m - p + 1)

  # scale factor
  scale_F <- c_corr * (m - p + 1) / (p * m)

  # threshold inverse scaling due to the scaled RD stated in Hardin&Rocke(2005)
  q_sHR <- stats::qf(p = 1 - quantile, df1 = df_F[1], df2 = df_F[2]) / scale_F

  flagged_outliers <- which(RD_obj$RD > q_sHR)

  result <- list(
    c = c_corr,
    m = m,
    df = df_F,
    scale = scale_F,
    threshold = q_sHR,
    flagged_outliers = flagged_outliers,
    call = call
  )

  class(result) <- if (SHASH) "SHASH_result" else "F_result"
  result
}
