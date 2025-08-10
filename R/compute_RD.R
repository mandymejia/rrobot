#' Compute Robust Distance and Covariance from a Subset
#'
#' Calculates the robust mean, covariance matrix, and optionally robust distances
#' using either:
#' - "auto" mode: automatically selects the best robust subset using covMcd
#' - "manual" mode: uses provided robust covariance matrix and subset indices
#'
#' @param x A numeric matrix or data frame of dimensions T × Q (observations × variables).
#' @param mode Character string; either "auto" (default) to compute MCD internally or "manual" to use user-supplied values.
#' @param cov_mcd Optional covariance matrix (Q × Q); required in "manual" mode.
#' @param ind_incld Optional vector of row indices used to compute the robust mean; required in "manual" mode.
#' @param dist Logical; if TRUE, compute squared robust Mahalanobis distances for all observations.
#'
#' @return A list with elements:
#' \describe{
#'   \item{ind_incld}{Vector of row indices used to compute the robust mean and covariance.}
#'   \item{ind_excld}{Vector of excluded row indices.}
#'   \item{h}{Number of included observations.}
#'   \item{xbar_star}{Robust mean vector (length Q).}
#'   \item{S_star}{Robust covariance matrix (Q × Q).}
#'   \item{invcov_sqrt}{Matrix square root of the inverse covariance matrix (Q × Q).}
#'   \item{RD}{Squared robust distances for all observations (length T), or NULL if dist = FALSE.}
#'   \item{call}{The matched function call.}
#' }
#'
#' @importFrom expm sqrtm
#' @importFrom robustbase covMcd
#' @importFrom stats cov
#' @export
compute_RD <- function(x, mode = c("auto", "manual"),
                    cov_mcd = NULL, ind_incld = NULL, dist = TRUE) {

  call <- match.call()

  mode <- match.arg(mode)
  x <- as.matrix(x)
  stopifnot(is.numeric(x))

  t <- nrow(x)
  Q <- ncol(x)

  if (mode == "auto") {
    cov_obj <- robustbase::covMcd(x)
    ind_incld <- cov_obj$best
    data_incld <- x[ind_incld, , drop = FALSE]
    cov_mcd <- cov(data_incld)
  } else {
    if (is.null(cov_mcd) || is.null(ind_incld)) {
      stop("In manual mode, both 'cov_mcd' and 'ind_incld' must be provided.")
    }
    stopifnot(all(ind_incld %in% seq_len(t)))
    data_incld <- x[ind_incld, , drop = FALSE]
  }

  h <- length(ind_incld)
  ind_excld <- setdiff(seq_len(t), ind_incld)

  if (is.null(dim(data_incld)) || ncol(data_incld) == 1) {
    data_incld <- matrix(data_incld, ncol = 1)
  }

  # Compute robust mean
  xbar_star <- colMeans(data_incld)

  # Defensive check for invertibility
  if (qr(cov_mcd)$rank < ncol(cov_mcd)) {
    warning("Covariance matrix may be rank-deficient; inverse may be unstable.")
  }

  invcov <- solve(cov_mcd)
  invcov_sqrt <- expm::sqrtm(invcov)
  if (is.null(invcov_sqrt)) stop("Covariance inversion failed")

  # Compute robust distances if requested
  if (dist) {
    xbar_star_mat <- matrix(xbar_star, nrow = t, ncol = Q, byrow = TRUE)
    temp <- (x - xbar_star_mat) %*% invcov_sqrt
    RD <- rowSums(temp^2)
  } else {
    RD <- NULL
  }

  result <- list(
    ind_incld = ind_incld,
    ind_excld = ind_excld,
    h = h,
    xbar_star = xbar_star,
    S_star = cov_mcd,
    invcov_sqrt = invcov_sqrt,
    RD = RD,
    call = call
  )

  class(result) <- "RD_result"
  result
}
