#' Compute Squared robust distance and covariance from a Subset
#'
#' Calculates the robust mean, covariance matrix, and optionally robust distances
#' using either:
#' - "auto" mode: automatically selects the best robust subset using covMcd
#' - "manual" mode: uses provided robust covariance matrix and subset indices
#'
#' @inheritParams x
#' @inheritParams mode
#' @inheritParams cov_mcd
#' @inheritParams ind_incld
#' @inheritParams dist
#'
#' @return A list with elements:
#' \describe{
#'   \item{ind_incld}{Vector of row indices used to compute the robust mean and covariance.}
#'   \item{ind_excld}{Vector of excluded row indices.}
#'   \item{h}{Number of included observations.}
#'   \item{xbar_star}{Robust mean vector (length p).}
#'   \item{S_star}{Robust covariance matrix (p × p).}
#'   \item{invcov_sqrt}{Matrix square root of the inverse covariance matrix (p × p).}
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
  mode <- match.arg(mode)
  x <- as.matrix(x)
  stopifnot(is.numeric(x))
  n <- nrow(x); p <- ncol(x)

  if (mode == "auto") {
    cov_obj   <- robustbase::covMcd(x)
    ind_incld <- cov_obj$best
    S_star    <- cov_obj$cov
    xbar_star <- as.numeric(cov_obj$center)
  } else {
    if (is.null(cov_mcd) || is.null(ind_incld))
      stop("In manual mode, both 'cov_mcd' and 'ind_incld' must be provided.")
    stopifnot(all(ind_incld %in% seq_len(n)))
    S_star    <- cov_mcd
    xbar_star <- colMeans(x[ind_incld, , drop = FALSE])
  }

  h <- length(ind_incld)
  ind_excld <- setdiff(seq_len(n), ind_incld)

  if (qr(S_star)$rank < p)
    warning("S* may be rank-deficient; distances could be unstable.")

  RD2_all <- if (dist) stats::mahalanobis(x, center = xbar_star, cov = S_star) else NULL

  result <- list(
    ind_incld  = ind_incld,
    ind_excld  = ind_excld,
    n = n, p = p, h = h,
    xbar_star  = xbar_star,
    S_star     = S_star,
    RD   = RD2_all,
    RD_incld  = if (!is.null(RD2_all)) RD2_all[ind_incld] else NULL,
    RD_excld  = if (!is.null(RD2_all)) RD2_all[ind_excld] else NULL,
    call = call
  )

  class(result) <- "RD_result"
  result
}
