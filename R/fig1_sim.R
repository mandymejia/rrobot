#' Simulate robust distances and F-fit thresholding
#'
#' Performs robust covariance estimation on the input data using MCD,
#' then fits a Hardin & Rocke F-distribution threshold, and counts
#' the number of observations exceeding this threshold.
#'
#' @param data A numeric matrix or data frame of dimensions \code{T x Q}.
#' @param alpha Significance level for the F-distribution threshold.
#'
#' @return A list containing:
#' \describe{
#'   \item{RD}{Vector of squared robust distances.}
#'   \item{ind_incld}{Indices of observations included in the robust covariance estimate.}
#'   \item{Y}{Original input data.}
#'   \item{scale}{Scale parameter of the fitted F-distribution.}
#'   \item{num_out}{Number of observations flagged as outliers.}
#'   \item{q_sHR}{The theoretical F-distribution quantile threshold.}
#' }
#' @export
fig1_sim <- function(data, alpha) {
  t <- nrow(data)
  Q <- ncol(data)

  rob_obj <- comp_RD(data)
  ffit <- Fit_F(Q, n = t, h = rob_obj$h, quantile = alpha)
  num_out <- length(which(rob_obj$RD >= ffit$threshold))

  list(
    RD = rob_obj$RD,
    ind_incld = rob_obj$ind_incld,
    Y = data,
    scale = ffit$scale,
    num_out = num_out,
    q_sHR = ffit$threshold
  )
}
