
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

#' Outlier Initialization
#'
#' Initializes outlier candidate weights using one of three methods:
#' isolation forest (isoplus), empirical rule, or the 50% distance rule.
#' Returns initial inlier weights and threshold cutoffs on the original scale.
#'
#' @param x Numeric vector.
#' @param method Character. One of \code{"isoplus"}, \code{"emprule"},
#'   or \code{"fifty"}.
#' @param thr Positive numeric scalar. Threshold for empirical rule
#'   (only used when \code{method = "emprule"}).
#' @param thr_iso Numeric scalar in \[0, 1\]. Isolation forest anomaly score
#'   threshold (only used when \code{method = "isoplus"}).
#' @param tail Character. One of \code{"both"}, \code{"upper"},
#'   or \code{"lower"}.
#' @param seed Integer or NULL. Random seed for isolation forest
#'   (only used when \code{method = "isoplus"}).
#' @param ntrees Integer. Number of trees for isolation forest
#'   (only used when \code{method = "isoplus"}).
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{weight_new}}{Logical vector. TRUE = inlier.}
#'     \item{\code{iso_scores}}{Numeric vector of isolation forest scores,
#'       or NA if not used.}
#'     \item{\code{cutoff_upper}}{Upper cutoff on original scale, or NA.}
#'     \item{\code{cutoff_lower}}{Lower cutoff on original scale, or NA.}
#'   }
#'
#' @importFrom isotree isolation.forest
#' @importFrom stats predict median
#' @export
outlier_init <- function(x,
                         method  = c("isoplus", "emprule", "fifty"),
                         thr     = 2.58,
                         thr_iso = 0.6,
                         tail    = c("both", "upper", "lower"),
                         seed    = NULL,
                         ntrees  = 200) {

  method <- match.arg(method)
  tail   <- match.arg(tail)

  x_med        <- median(x, na.rm = TRUE)
  iso_scores   <- rep(NA_real_, length(x))
  cutoff_upper <- NA_real_
  cutoff_lower <- NA_real_

  # -----------------------------------------------------------
  # Method 1: isoplus — Isolation Forest Initialization
  # -----------------------------------------------------------
  if (method == "isoplus") {

    iso_model  <- isotree::isolation.forest(
      data.frame(x),
      ntrees   = ntrees,
      seed     = if (!is.null(seed)) seed else 1
    )
    iso_scores <- predict(iso_model, data.frame(x), type = "score")
    idx_out0   <- which(iso_scores >= thr_iso)

    if (length(idx_out0) == 0) {
      weight_new <- rep(TRUE, length(x))
    } else {
      upp <- idx_out0[x[idx_out0] >= x_med]
      low <- idx_out0[x[idx_out0] <= x_med]

      if (tail == "upper") {
        cutoff_upper <- if (length(upp) == 0) Inf else min(x[upp], na.rm = TRUE)
        weight_new   <- x < cutoff_upper
      } else if (tail == "lower") {
        cutoff_lower <- if (length(low) == 0) -Inf else max(x[low], na.rm = TRUE)
        weight_new   <- x > cutoff_lower
      } else {
        cutoff_upper <- if (length(upp) == 0) Inf  else min(x[upp], na.rm = TRUE)
        cutoff_lower <- if (length(low) == 0) -Inf else max(x[low], na.rm = TRUE)
        weight_new   <- (x > cutoff_lower) & (x < cutoff_upper)
      }
    }

    # -----------------------------------------------------------
    # Method 2: emprule — Robust Empirical Rule Initialization
    # -----------------------------------------------------------
  } else if (method == "emprule") {

    flagged    <- emprule_rob(x, thr = thr, tail = tail)
    weight_new <- !flagged

    mad_val <- 1.4826 * median(abs(x - x_med), na.rm = TRUE)

    if (tail %in% c("both", "upper")) {
      cutoff_upper <- x_med + thr * mad_val
    }
    if (tail %in% c("both", "lower")) {
      cutoff_lower <- x_med - thr * mad_val
    }

    # -----------------------------------------------------------
    # Method 3: fifty — Top 50% by distance from median
    # -----------------------------------------------------------
  } else if (method == "fifty") {

    dist_from_med <- abs(x - x_med)

    if (tail == "upper") {
      candidates <- which(x > x_med)
      n_remove   <- floor(length(candidates) * 0.5)
      ranked     <- candidates[order(dist_from_med[candidates], decreasing = TRUE)]
      remove_idx <- ranked[seq_len(n_remove)]
      weight_new <- rep(TRUE, length(x))
      weight_new[remove_idx] <- FALSE
      kept         <- setdiff(candidates, remove_idx)
      cutoff_upper <- if (length(kept) > 0) max(x[kept], na.rm = TRUE) else x_med

    } else if (tail == "lower") {
      candidates <- which(x < x_med)
      n_remove   <- floor(length(candidates) * 0.5)
      ranked     <- candidates[order(dist_from_med[candidates], decreasing = TRUE)]
      remove_idx <- ranked[seq_len(n_remove)]
      weight_new <- rep(TRUE, length(x))
      weight_new[remove_idx] <- FALSE
      kept         <- setdiff(candidates, remove_idx)
      cutoff_lower <- if (length(kept) > 0) min(x[kept], na.rm = TRUE) else x_med

    } else {
      n_total    <- length(x)
      n_remove   <- floor(n_total * 0.5)
      ranked     <- order(dist_from_med, decreasing = TRUE)
      remove_idx <- ranked[seq_len(n_remove)]
      weight_new <- rep(TRUE, n_total)
      weight_new[remove_idx] <- FALSE
      kept         <- which(weight_new)
      cutoff_upper <- max(x[kept], na.rm = TRUE)
      cutoff_lower <- min(x[kept], na.rm = TRUE)
    }
  }

  list(
    weight_new   = weight_new,
    iso_scores   = iso_scores,
    cutoff_upper = cutoff_upper,
    cutoff_lower = cutoff_lower
  )
}


# =============================================================
# SHASH_out
# =============================================================

#' SHASH-based Outlier Detection (Extended)
#'
#' Detects univariate outliers using an iterative SHASH fitting process.
#' The initialization method is controlled via \code{method_init}.
#'
#' @param x Numeric vector. May contain \code{NA} values.
#' @param thr0 Positive numeric scalar. Threshold passed to
#'   \code{outlier_init} for empirical rule initialization (default: 2.58).
#' @param thr1 Positive numeric scalar. Threshold for iterative convergence
#'   (default: 2.58).
#' @param thr Positive numeric scalar. Final outlier detection threshold
#'   (default: 4).
#' @param tail Character. One of \code{"both"}, \code{"upper"},
#'   or \code{"lower"}.
#' @param method_init Character. Initialization method. One of
#'   \code{"isoplus"} (default), \code{"emprule"}, or \code{"fifty"}.
#' @param thr_iso Numeric scalar in \[0, 1\]. Isolation forest threshold
#'   (default: 0.6). Only used when \code{method_init = "isoplus"}.
#' @param iso_seed Integer or NULL. Seed for isolation forest. Only used
#'   when \code{method_init = "isoplus"}.
#' @param maxit Positive integer. Maximum iterations (default: 100).
#' @param weight_init Optional logical vector. Bypasses initialization
#'   if supplied.
#'
#' @return A list of class \code{"SHASH_out"} with elements:
#'   \describe{
#'     \item{\code{out_idx}}{Integer vector of outlier indices.}
#'     \item{\code{x_norm}}{Numeric vector of SHASH-normalised scores.}
#'     \item{\code{SHASH_coef}}{Named list of fitted SHASH parameters.}
#'     \item{\code{isotree_scores}}{Isolation forest scores, or NA.}
#'     \item{\code{initial_weights}}{Logical vector of initialization weights.}
#'     \item{\code{init_cutoff_upper}}{Upper initialization cutoff on original scale.}
#'     \item{\code{init_cutoff_lower}}{Lower initialization cutoff on original scale.}
#'     \item{\code{iter_thr1_upper}}{thr1 back-transformed per iteration (upper).}
#'     \item{\code{iter_thr1_lower}}{thr1 back-transformed per iteration (lower).}
#'     \item{\code{final_thr_upper}}{Final thr back-transformed (upper).}
#'     \item{\code{final_thr_lower}}{Final thr back-transformed (lower).}
#'     \item{\code{indx_iters}}{Integer matrix of outlier flags per iteration.}
#'     \item{\code{norm_iters}}{Numeric matrix of normalised scores per iteration.}
#'     \item{\code{last_iter}}{Number of iterations completed.}
#'     \item{\code{converged}}{Logical. TRUE if converged before maxit.}
#'     \item{\code{params}}{List of input parameters.}
#'   }
#'
#' @examples
#' # --- Example 1: SHASH-i (isoplus initialization) ---
#' set.seed(42)
#' x <- rnorm(200, mean = 10, sd = 2)
#' x[c(17, 77)] <- x[c(17, 77)] + 10
#'
#' res <- SHASH_out(x, thr0 = 2.58, thr1 = 2.58, thr = 3,
#'                  tail = "both", method_init = "isoplus")
#' res$out_idx
#' res$converged
#'
#' # --- Example 2: SHASH-z (emprule initialization) ---
#' res_z <- SHASH_out(x, thr0 = 2.58, thr1 = 2.58, thr = 3,
#'                    tail = "both", method_init = "emprule")
#' res_z$out_idx
#'
#' # --- Example 3: SHASH-50 (fifty initialization) ---
#' res_50 <- SHASH_out(x, thr0 = 2.58, thr1 = 2.58, thr = 3,
#'                     tail = "both", method_init = "fifty")
#' res_50$out_idx
#'
#' # --- Example 4: Real benchmark data (HBK) ---
#' data("hbk", package = "robustbase")
#' res_hbk <- SHASH_out(hbk$X1, thr0 = 2.58, thr1 = 2.58, thr = 3,
#'                      tail = "both", method_init = "isoplus")
#' res_hbk$out_idx
#' res_hbk$converged
#'
#' @importFrom gamlss gamlssML coefAll
#' @export
SHASH_out <- function(x,
                      thr0        = 2.58,
                      thr1        = 2.58,
                      thr         = 4,
                      tail        = c("both", "upper", "lower"),
                      method_init = c("isoplus", "emprule", "fifty"),
                      thr_iso     = 0.6,
                      iso_seed    = NULL,
                      maxit       = 100,
                      weight_init = NULL) {

  stopifnot("Input 'x' must be a numeric vector, not a matrix" = !is.matrix(x))
  tail        <- match.arg(tail)
  method_init <- match.arg(method_init)

  if (method_init == "isoplus" && (thr_iso < 0 || thr_iso > 1)) {
    warning("`thr_iso` should be between 0 and 1 when using isolation forest.")
  }
  if (!is.numeric(thr_iso) || length(thr_iso) != 1 || thr_iso < 0 || thr_iso > 1) {
    stop("'thr_iso' must be a numeric value between 0 and 1.")
  }
  if (!is.null(weight_init) && length(weight_init) != length(x)) {
    stop("'weight_init' must have the same length as 'x'.")
  }

  params <- list(
    orig_values  = x,
    thr0         = thr0,
    thr1         = thr1,
    thr          = thr,
    tail         = tail,
    method_init  = method_init,
    thr_iso      = thr_iso,
    iso_seed     = iso_seed,
    maxit        = maxit,
    weight_init  = weight_init
  )

  # helpers for tail logic
  inlier_by_tail <- function(z, thr, tail) {
    if (tail == "upper")      z <= thr
    else if (tail == "lower") z >= -thr
    else                      abs(z) <= thr
  }

  outlier_idx_by_tail <- function(z, thr, tail) {
    if (tail == "upper")      which(z > thr)
    else if (tail == "lower") which(z < -thr)
    else                      which(abs(z) > thr)
  }

  na_locs <- is.na(x)
  x_clean <- x[!na_locs]
  n_clean <- length(x_clean)

  # ---------------------------------------------------------
  # Initialization
  # ---------------------------------------------------------

  iso_scores        <- rep(NA_real_, n_clean)
  init_cutoff_upper <- NA_real_
  init_cutoff_lower <- NA_real_

  if (!is.null(weight_init)) {
    weight_new <- as.logical(weight_init[!na_locs])

  } else {
    init_result       <- outlier_init(
      x_clean,
      method  = method_init,
      thr     = thr0,
      thr_iso = thr_iso,
      tail    = tail,
      seed    = iso_seed
    )
    weight_new        <- init_result$weight_new
    iso_scores        <- init_result$iso_scores
    init_cutoff_upper <- init_result$cutoff_upper
    init_cutoff_lower <- init_result$cutoff_lower
  }

  initial_weights_clean <- weight_new

  # ---------------------------------------------------------
  # Iteration loop
  # ---------------------------------------------------------

  norm_iters_clean <- matrix(NA_real_,    nrow = n_clean, ncol = maxit)
  indx_iters_clean <- matrix(NA_integer_, nrow = n_clean, ncol = maxit)
  iter_thr1_upper  <- rep(NA_real_, maxit)
  iter_thr1_lower  <- rep(NA_real_, maxit)

  iter    <- 0
  success <- FALSE
  est     <- NULL

  repeat {
    iter       <- iter + 1
    weight_old <- weight_new

    mod <- gamlss::gamlssML(
      x_clean ~ 1, family = "SHASHo2", maxit = 1e4,
      weights = as.numeric(weight_new)
    )
    est <- gamlss::coefAll(mod)

    x_norm_clean <- SHASH_to_normal(
      x_clean,
      mu = est$mu, sigma = est$sigma,
      nu = est$nu, tau   = est$tau
    )
    norm_iters_clean[, iter] <- x_norm_clean

    if (tail %in% c("both", "upper")) {
      iter_thr1_upper[iter] <- tryCatch(
        normal_to_SHASH(thr1,  mu = est$mu, sigma = est$sigma,
                        nu = est$nu, tau   = est$tau),
        error = function(e) NA_real_
      )
    }
    if (tail %in% c("both", "lower")) {
      iter_thr1_lower[iter] <- tryCatch(
        normal_to_SHASH(-thr1, mu = est$mu, sigma = est$sigma,
                        nu = est$nu, tau   = est$tau),
        error = function(e) NA_real_
      )
    }

    weight_new               <- inlier_by_tail(x_norm_clean, thr1, tail)
    indx_iters_clean[, iter] <- as.integer(!weight_new)

    if (isTRUE(all.equal(weight_old, weight_new))) {
      success <- TRUE
      break
    }
    if (iter >= maxit) break
  }

  norm_iters_clean <- norm_iters_clean[, seq_len(iter), drop = FALSE]
  indx_iters_clean <- indx_iters_clean[, seq_len(iter), drop = FALSE]
  iter_thr1_upper  <- iter_thr1_upper[seq_len(iter)]
  iter_thr1_lower  <- iter_thr1_lower[seq_len(iter)]

  # ---------------------------------------------------------
  # Final threshold back-transformed
  # ---------------------------------------------------------

  final_thr_upper <- if (tail %in% c("both", "upper")) {
    tryCatch(
      normal_to_SHASH(thr,  mu = est$mu, sigma = est$sigma,
                      nu = est$nu, tau   = est$tau),
      error = function(e) NA_real_
    )
  } else NA_real_

  final_thr_lower <- if (tail %in% c("both", "lower")) {
    tryCatch(
      normal_to_SHASH(-thr, mu = est$mu, sigma = est$sigma,
                      nu = est$nu, tau   = est$tau),
      error = function(e) NA_real_
    )
  } else NA_real_

  # ---------------------------------------------------------
  # Pad back to original length
  # ---------------------------------------------------------

  n_orig  <- length(x)
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
  out_flag_full[!na_locs][outlier_idx_by_tail(x_norm_clean, thr, tail)] <- TRUE
  final_out_idx_full <- which(out_flag_full)

  out <- list(
    out_idx           = final_out_idx_full,
    x_norm            = x_norm_full,
    SHASH_coef        = est[c("mu", "sigma", "nu", "tau")],
    isotree_scores    = isotree_scores_full,
    initial_weights   = initial_weights_full,
    init_cutoff_upper = init_cutoff_upper,
    init_cutoff_lower = init_cutoff_lower,
    iter_thr1_upper   = iter_thr1_upper,
    iter_thr1_lower   = iter_thr1_lower,
    final_thr_upper   = final_thr_upper,
    final_thr_lower   = final_thr_lower,
    indx_iters        = indx_iters_full,
    norm_iters        = norm_iters_full,
    last_iter         = iter,
    converged         = success,
    params            = params
  )
  class(out) <- "SHASH_out"
  out
}

# -----------------------------------------------------------------------------
# Robust Empirical Rule
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

