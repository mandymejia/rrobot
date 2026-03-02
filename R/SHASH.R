
#' SHASH Data Transformation
#'
#' @description
#' These two functions form a matched pair for transforming data between the
#' SHASH (Sinh-Arcsinh) distribution and the standard normal distribution.
#' `SHASH_to_normal()` maps SHASH-distributed observations onto an
#' approximately normal scale; `normal_to_SHASH()` is the inverse.
#'
#' @examples
#' set.seed(42)
#' x <- rnorm(200)
#' x[c(17, 77)] <- x[c(17, 77)] + 5
#'
#' mu <- 0; sigma <- 0; nu <- 0; tau <- 0
#'
#' z <- SHASH_to_normal(x, mu = mu, sigma = sigma, nu = nu, tau = tau)
#' x_recovered <- normal_to_SHASH(z, mu = mu, sigma = sigma, nu = nu, tau = tau)
#' all.equal(x, x_recovered)
#'
#' @name SHASH_transform
NULL


#' @describeIn SHASH_transform Transforms SHASH-distributed data to
#'   approximately normal data.
#'
#' @param x Numeric vector of values to transform.
#' @param mu Numeric scalar. Location parameter controlling the mean of the
#'   SHASH distribution.
#' @param sigma Numeric scalar. Spread parameter on the log scale. The function
#'   applies \code{exp(sigma)} internally, so pass the raw coefficient as
#'   returned by \code{gamlssML()}. Pass \code{sigma = 0} to get unit spread
#'   since \code{exp(0) = 1}.
#' @param nu Numeric scalar. Skewness parameter. A value of 0 gives a
#'   symmetric distribution.
#' @param tau Numeric scalar. Tail-weight parameter on the log scale. Pass
#'   \code{tau = 0} for normal-like tails since \code{exp(0) = 1}.
#'
#' @return A numeric vector of transformed values, the same length as \code{x}.
#' @export
SHASH_to_normal <- function(x, mu, sigma, nu, tau){
  stopifnot(
    "x must be a numeric vector"           = is.numeric(x),
    "mu must be a single numeric value"    = is.numeric(mu)    && length(mu)    == 1,
    "sigma must be a single numeric value" = is.numeric(sigma) && length(sigma) == 1,
    "nu must be a single numeric value"    = is.numeric(nu)    && length(nu)    == 1,
    "tau must be a single numeric value"   = is.numeric(tau)   && length(tau)   == 1
  )
  sigma <- exp(sigma)
  tau   <- exp(tau)
  sinh((tau * asinh((x - mu) / (sigma * tau))) - nu)
}


#' @describeIn SHASH_transform Transforms standard normal data back to the
#'   SHASH-distributed scale.
#' @export
normal_to_SHASH <- function(x, mu, sigma, nu, tau){
  stopifnot(
    "x must be a numeric vector"           = is.numeric(x),
    "mu must be a single numeric value"    = is.numeric(mu)    && length(mu)    == 1,
    "sigma must be a single numeric value" = is.numeric(sigma) && length(sigma) == 1,
    "nu must be a single numeric value"    = is.numeric(nu)    && length(nu)    == 1,
    "tau must be a single numeric value"   = is.numeric(tau)   && length(tau)   == 1
  )
  sigma <- exp(sigma)
  tau   <- exp(tau)
  (sigma * tau * sinh((asinh(x) + nu) / tau)) + mu
}

#' SHASH-based Outlier Detection (Extended)
#'
#' Detects univariate outliers using an iterative SHASH fitting process with
#' optional pre-flagging strategies. A SHASH (Sinh-Arcsinh) distribution is
#' fitted to the data iteratively, each time excluding candidate outliers from
#' the fit, until the set of flagged observations converges or \code{maxit}
#' is reached.
#'
#' @param x Numeric vector. May contain \code{NA} values; they are excluded
#'   from fitting and propagated as \code{NA} in all output vectors.
#' @param thr0 Positive numeric scalar. Threshold for initial outlier
#'   pre-flagging when \code{use_iso = FALSE} (default: 2.58).
#' @param thr1 Positive numeric scalar. Threshold used to classify observations
#'   as inliers during iterative convergence (default: 2.58).
#' @param thr Positive numeric scalar. Final threshold applied to the converged
#'   SHASH-normalised scores to declare outliers in the returned output
#'   (default: 4).
#' @param tail Character string specifying which tail(s) to check for outliers.
#'   Must be one of \code{"both"} (default), \code{"upper"}, or
#'   \code{"lower"}.
#'   \itemize{
#'     \item \code{"upper"}: detect upper-tail outliers only.
#'     \item \code{"lower"}: detect lower-tail outliers only.
#'     \item \code{"both"}: detect two-sided outliers.
#'   }
#' @param use_iso Logical. If \code{TRUE} (default), uses an isolation forest
#'   (via \code{isotree}) to pre-screen candidate outliers before the iterative
#'   fitting loop begins.
#' @param thr_iso Numeric scalar in \[0, 1\]. Isolation forest anomaly score
#'   threshold above which observations are treated as candidate outliers
#'   during pre-screening (default: 0.6). Only used when \code{use_iso = TRUE}.
#' @param maxit Positive integer. Maximum number of fitting iterations before
#'   the algorithm stops regardless of convergence (default: 100).
#' @param weight_init Optional logical vector of length \code{length(x)}.
#'   If supplied, these weights initialise the iterative fit directly,
#'   bypassing both the isolation forest and empirical-rule pre-screening.
#'   \code{TRUE} means the observation is treated as an inlier in the first
#'   iteration.
#'
#' @return A list of class \code{"SHASH_out"} with the following elements:
#'   \describe{
#'     \item{\code{out_idx}}{Integer vector. Indices of observations in
#'       \code{x} that were flagged as outliers at the final threshold
#'       \code{thr}.}
#'     \item{\code{x_norm}}{Numeric vector. SHASH-normalised scores for every
#'       observation (same length as \code{x}; \code{NA} where \code{x} was
#'       \code{NA}).}
#'     \item{\code{SHASH_coef}}{Named list with elements \code{mu},
#'       \code{sigma}, \code{nu}, and \code{tau}: the fitted SHASH parameter
#'       estimates from the final iteration (sigma and tau are on the log
#'       scale, as returned by \code{gamlssML}).}
#'     \item{\code{isotree_scores}}{Numeric vector of isolation forest anomaly
#'       scores (same length as \code{x}). \code{NA} when \code{use_iso =
#'       FALSE} or \code{weight_init} was supplied.}
#'     \item{\code{initial_weights}}{Logical vector. Inlier weights used for
#'       the very first fitting iteration (same length as \code{x}).}
#'     \item{\code{indx_iters}}{Integer matrix of dimensions
#'       \code{length(x)} × \code{last_iter}. Each column records which
#'       observations were flagged as outliers (value 1) during that iteration.}
#'     \item{\code{norm_iters}}{Numeric matrix of dimensions
#'       \code{length(x)} × \code{last_iter}. Each column records the
#'       SHASH-normalised scores from that iteration.}
#'     \item{\code{last_iter}}{Integer. The number of iterations completed
#'       before convergence or hitting \code{maxit}.}
#'     \item{\code{converged}}{Logical. \code{TRUE} if the inlier weight
#'       vector stabilised before reaching \code{maxit}.}
#'     \item{\code{params}}{List. A record of all input parameters, stored
#'       for reproducibility.}
#'   }
#'
#' @importFrom isotree isolation.forest
#' @importFrom stats predict
#' @importFrom gamlss gamlssML coefAll
#'
#' @examples
#' # --- Example 1: Synthetic data with known injected outliers ---------------
#' # Using rnorm lets us inject outliers at known positions so we can verify
#' # the function finds exactly what we planted.
#' set.seed(42)
#' x <- rnorm(200, mean = 10, sd = 2)
#'
#' # Shift a handful of observations far into the upper tail
#' outlier_positions <- c(17, 77, seq(190, 200))
#' x[outlier_positions] <- x[outlier_positions] + 10
#'
#' result_sim <- SHASH_out(
#'   x,
#'   thr0    = 2.58,
#'   thr1    = 2.58,
#'   thr     = 4,
#'   tail    = "both",
#'   use_iso = FALSE   # skip isolation forest to keep the example fast
#' )
#'
#' result_sim$out_idx    # should recover positions near outlier_positions
#' result_sim$converged  # did the iterative fit stabilise?
#'
#' # --- Example 2: Real benchmark data (Hawkins-Bradu-Kass) ------------------
#' # hbk is a classic outlier detection benchmark shipped with robustbase,
#' # which this package already imports, so it is always available.
#' data("hbk", package = "robustbase")
#'
#' result_hbk <- SHASH_out(
#'   hbk$X1,
#'   thr0    = 2.58,
#'   thr1    = 2.58,
#'   thr     = 4,
#'   tail    = "both",
#'   use_iso = FALSE
#' )
#'
#' result_hbk$out_idx   # flagged observations in the X1 column
#' result_hbk$SHASH_coef  # fitted SHASH parameters; sigma and tau are log-scale
#'
#' # Which positions were flagged as outliers?
#' result_hbk$out_idx
#'
#' # Did the algorithm converge before hitting maxit?
#' result_hbk$converged
#'
#' # How many iterations did it take?
#' result_hbk$last_iter
#'
#' @export
SHASH_out <- function(x,
                      thr0         = 2.58,
                      thr1         = 2.58,
                      thr          = 4,
                      tail         = c("both", "upper", "lower"),
                      use_iso      = TRUE,
                      thr_iso      = 0.6,
                      maxit        = 100,
                      weight_init  = NULL) {

  stopifnot("Input 'x' must be a numeric vector, not a matrix" = !is.matrix(x))
  tail <- match.arg(tail)

  if (use_iso && (thr_iso < 0 || thr_iso > 1)) {
    warning("`thr_iso` should be between 0 and 1 when using isolation forest.")
  }
  if (!is.numeric(thr_iso) || length(thr_iso) != 1 || thr_iso < 0 || thr_iso > 1) {
    stop("Threshold Isolation Forest 'thr_iso' must be a positive numeric value between 0 and 1.")
  }

  # Validate weight_init length if supplied, to catch mismatches early
  # rather than silently recycling or truncating during subsetting.
  if (!is.null(weight_init) && length(weight_init) != length(x)) {
    stop("'weight_init' must have the same length as 'x'.")
  }

  params <- list(
    orig_values  = x,
    thr0         = thr0,
    thr1         = thr1,
    thr          = thr,
    tail         = tail,
    use_iso      = use_iso,
    thr_iso      = thr_iso,
    maxit        = maxit,
    weight_init  = weight_init
  )

  # helpers for tail logic
  inlier_by_tail <- function(z, thr, tail) {
    if (tail == "upper") {
      z <= thr
    } else if (tail == "lower") {
      z >= -thr
    } else {
      abs(z) <= thr
    }
  }
  outlier_idx_by_tail <- function(z, thr, tail) {
    if (tail == "upper") {
      which(z > thr)
    } else if (tail == "lower") {
      which(z < -thr)
    } else {
      which(abs(z) > thr)
    }
  }

  na_locs   <- is.na(x)
  x_clean   <- x[!na_locs]
  n_clean   <- length(x_clean)
  x_med     <- median(x_clean)

  if (!is.null(weight_init)) {
    weight_new <- as.logical(weight_init[!na_locs])
    iso_scores <- rep(NA_real_, n_clean)
  } else if (use_iso) {
    iso_model  <- isotree::isolation.forest(data.frame(x_clean), ntrees = 200)
    iso_scores <- predict(iso_model, data.frame(x_clean), type = "score")
    idx_out0   <- which(iso_scores >= thr_iso)
    if (length(idx_out0) == 0) {
      weight_new <- rep(TRUE, n_clean)
    } else {
      upp <- idx_out0[x_clean[idx_out0] >= x_med]
      low <- idx_out0[x_clean[idx_out0] <= x_med]
      if (tail == "upper") {
        weight_new <- if (length(upp) == 0) rep(TRUE, n_clean) else x_clean < min(x_clean[upp], na.rm = TRUE)
      } else if (tail == "lower") {
        weight_new <- if (length(low) == 0) rep(TRUE, n_clean) else x_clean > max(x_clean[low], na.rm = TRUE)
      } else {
        cutoff_upper <- if (length(upp) == 0) Inf  else min(x_clean[upp], na.rm = TRUE)
        cutoff_lower <- if (length(low) == 0) -Inf else max(x_clean[low], na.rm = TRUE)
        weight_new <- (x_clean > cutoff_lower) & (x_clean < cutoff_upper)
      }
    }
  } else {
    iso_scores <- rep(NA_real_, n_clean)
    W <- tryCatch(
      1 - emprule_rob(x_clean, thr = thr0, tail = tail),
      error = function(e) {
        warning("Empirical rule failed. Defaulting to all TRUE weights.")
        rep(TRUE, n_clean)
      }
    )
    weight_new <- as.logical(W)
  }

  initial_weights_clean <- weight_new

  norm_iters_clean <- matrix(NA_real_,    nrow = n_clean, ncol = maxit)
  indx_iters_clean <- matrix(NA_integer_, nrow = n_clean, ncol = maxit)
  iter <- 0
  success <- FALSE

  repeat {
    iter <- iter + 1
    weight_old <- weight_new
    mod <- gamlss::gamlssML(
      x_clean ~ 1, family = "SHASHo2", maxit = 1e4,
      weights = as.numeric(weight_new)
    )
    est <- gamlss::coefAll(mod)

    x_norm_clean <- SHASH_to_normal(
      x_clean, mu = est$mu, sigma = est$sigma,
      nu = est$nu, tau = est$tau
    )
    norm_iters_clean[, iter] <- x_norm_clean

    weight_new <- inlier_by_tail(x_norm_clean, thr1, tail)
    indx_iters_clean[, iter] <- as.integer(!weight_new)

    if (isTRUE(all.equal(weight_old, weight_new))) {
      success <- TRUE
      break
    }
    if (iter >= maxit) break
  }

  norm_iters_clean <- norm_iters_clean[, seq_len(iter), drop = FALSE]
  indx_iters_clean <- indx_iters_clean[, seq_len(iter), drop = FALSE]

  final_idx_clean <- outlier_idx_by_tail(x_norm_clean, thr, tail)

  n_orig <- length(x)
  pad_vec <- function(v) {
    v_full <- rep(NA, n_orig)
    v_full[!na_locs] <- v
    v_full
  }
  x_norm_full          <- pad_vec(x_norm_clean)
  initial_weights_full <- pad_vec(initial_weights_clean)
  isotree_scores_full  <- pad_vec(iso_scores)
  norm_iters_full      <- matrix(NA_real_,    nrow = n_orig, ncol = iter)
  norm_iters_full[!na_locs, ] <- norm_iters_clean
  indx_iters_full      <- matrix(NA_integer_, nrow = n_orig, ncol = iter)
  indx_iters_full[!na_locs, ] <- indx_iters_clean

  out_flag_full <- rep(FALSE, n_orig)
  out_flag_full[!na_locs][final_idx_clean] <- TRUE
  final_out_idx_full <- which(out_flag_full)

  out <- list(
    out_idx         = final_out_idx_full,
    x_norm          = x_norm_full,
    SHASH_coef      = est[c("mu","sigma","nu","tau")],
    isotree_scores  = isotree_scores_full,
    initial_weights = initial_weights_full,
    indx_iters      = indx_iters_full,
    norm_iters      = norm_iters_full,
    last_iter       = iter,
    converged       = success,
    params          = params
  )
  class(out) <- "SHASH_out"
  out
}


# -----------------------------------------------------------------------------
# emprule_rob (internal helper — not exported
# -----------------------------------------------------------------------------

#' Robust Empirical Rule Outlier Detection
#'
#' Detects outliers using the median ± \code{thr} × MAD rule, where MAD is
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
#' @keywords internal
emprule_rob <- function(x, thr = 4, tail = c("both", "upper", "lower")) {

  tail <- match.arg(tail)

  x_med  <- median(x, na.rm = TRUE)
  mad_val <- 1.4826 * median(abs(x - x_med), na.rm = TRUE)

  # If MAD is zero (e.g. more than half the data are identical), the rule
  # cannot distinguish outliers — return a safe all-FALSE vector.
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

