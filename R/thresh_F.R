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
#' @param Q Integer. The number of variables (dimension of the data).
#' @param n Integer. The total sample size.
#' @param h Integer. The number of observations retained in the MCD subset.
#' @param quantile Numeric in (0,1) specifying the upper quantile for thresholding.
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
thresh_F <- function(Q, n, h, quantile) {
  call <- match.call()

  # Step 1: consistency correction factor
  c <- pchisq(q = qchisq(p = h / n, df = Q), df = Q + 2) / (h / n)

  # Step 2: asymptotic approximation for degrees of freedom
  alpha <- (n - h) / n
  q_alpha <- qchisq(p = 1 - alpha, df = Q)

  c_alpha <- (1 - alpha) / pchisq(q = q_alpha, df = Q + 2)
  c2 <- -0.5 * pchisq(q = q_alpha, df = Q + 2)
  c3 <- -0.5 * pchisq(q = q_alpha, df = Q + 4)
  c4 <- 3 * c3

  b1 <- c_alpha * (c3 - c4) / (1 - alpha)
  b2 <- 0.5 + c_alpha * (c3 - (q_alpha / Q) * (c2 + (1 - alpha) / 2)) / (1 - alpha)

  v1 <- (1 - alpha) * b1^2 * (alpha * (c_alpha * q_alpha / Q - 1)^2 - 1) -
    2 * c3 * c_alpha^2 * (3 * (b1 - Q * b2)^2 + (Q + 2) * b2 * (2 * b1 - Q * b2))

  v2 <- n * (b1 * (b1 - Q * b2) * (1 - alpha))^2 * c_alpha^2
  v <- v1 / v2

  m <- 2 / (c_alpha^2 * v)

  # Step 3: finite-sample correction
  m <- m * exp(0.725 - 0.00663 * Q - 0.078 * log(n))

  df <- c(Q, m - Q + 1)

  # scale factor
  scale <- c * (m - Q + 1) / (Q * m)

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

  class(result) <- "F_result"
  result
}
