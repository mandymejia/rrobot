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
#' @inheritParams quantile
#' @param SHASH Boolean. If running SHASH variant.
#' @inheritParams verbose
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
#' @importFrom stats pchisq qchisq qf
#' @export
thresh_F <- function(p, n, h, quantile, SHASH = FALSE, verbose = FALSE) {
  if (verbose) message("Running ", if (SHASH) "SHASH" else "F", " method: F-distribution threshold...")
  call <- match.call()

  # Step 1: consistency correction factor
  c <- pchisq(q = qchisq(p = h / n, df = p), df = p + 2) / (h / n)

  # Step 2: asymptotic approximation for degrees of freedom
  alpha <- (n - h) / n
  q_alpha <- qchisq(p = 1 - alpha, df = p)

  c_alpha <- (1 - alpha) / pchisq(q = q_alpha, df = p + 2)
  c2 <- -0.5 * pchisq(q = q_alpha, df = p + 2)
  c3 <- -0.5 * pchisq(q = q_alpha, df = p + 4)
  c4 <- 3 * c3

  b1 <- c_alpha * (c3 - c4) / (1 - alpha)
  b2 <- 0.5 + c_alpha * (c3 - (q_alpha / p) * (c2 + (1 - alpha) / 2)) / (1 - alpha)

  v1 <- (1 - alpha) * b1^2 * (alpha * (c_alpha * q_alpha / p - 1)^2 - 1) -
    2 * c3 * c_alpha^2 * (3 * (b1 - p * b2)^2 + (p + 2) * b2 * (2 * b1 - p * b2))

  v2 <- n * (b1 * (b1 - p * b2) * (1 - alpha))^2 * c_alpha^2
  v <- v1 / v2

  m <- 2 / (c_alpha^2 * v)

  # Step 3: finite-sample correction
  m <- m * exp(0.725 - 0.00663 * p - 0.078 * log(n))

  df <- c(p, m - p + 1)

  # scale factor
  scale <- c * (m - p + 1) / (p * m)

  # threshold inverse scaling due to the scaled RD stated in Hardin&Rocke(2005)
  q_sHR <- qf(p = 1 - quantile, df1 = df[1], df2 = df[2]) / scale

  result <- list(
    c = c,
    m = m,
    df = df,
    scale = scale,
    threshold = q_sHR,
    call = call
  )

  if (SHASH) {
    class(result) <- "SHASH_result"
  } else {
    class(result) <- "F_result"
  }

  result
}
