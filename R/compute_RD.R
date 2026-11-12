#' Compute Squared robust distance and covariance from a Subset
#'
#' Calculates the robust mean, covariance matrix, and optionally robust distances
#' using either:
#' - "auto" mode: automatically selects the best robust subset using covMcd
#' - "manual" mode: uses provided robust covariance matrix and subset indices
#'
#' @inheritParams params_doc
#'
#' @return A list of class \code{"RD_result"} with elements:
#' \describe{
#'   \item{ind_incld}{Vector of row indices used to compute the robust mean and covariance.}
#'   \item{ind_excld}{Vector of excluded row indices.}
#'   \item{n}{Number of observations (rows of \code{x}).}
#'   \item{p}{Number of variables (columns of \code{x}).}
#'   \item{h}{Number of included observations.}
#'   \item{xbar_star}{Robust mean vector (length p).}
#'   \item{S_star}{Robust covariance matrix (p x p).}
#'   \item{RD}{Squared robust distances for all observations (length T), or NULL if dist = FALSE.}
#'   \item{RD_incld}{Squared robust distances of the included observations, or NULL if dist = FALSE.}
#'   \item{RD_excld}{Squared robust distances of the excluded observations, or NULL if dist = FALSE.}
#'   \item{call}{The matched function call.}
#' }
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200), ncol = 2)
#' x[1:3, ] <- x[1:3, ] + 6
#' rd <- compute_RD(x)
#' which(rd$RD > stats::qchisq(0.99, df = 2))
#'
#' @importFrom robustbase covMcd
#' @export
compute_RD <- function(x, mode = c("auto", "manual"),
                       cov_mcd = NULL, ind_incld = NULL, dist = TRUE) {
  call <- match.call()
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
    stopifnot("'cov_mcd' must be a p x p matrix" =
                is.matrix(cov_mcd) && all(dim(cov_mcd) == p))
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
