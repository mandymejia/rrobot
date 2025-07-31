#' SHASH Data Transformation
#'
#' Transforms data between the SHASH (Sinh-Arcsinh) distribution and the standard normal distribution.
#'
#' @param x Numeric vector of data to transform.
#' @param mu Parameter controlling location (mean) of the SHASH distribution.
#' @param sigma Parameter controlling spread (variance); must be positive (on log scale).
#' @param nu Parameter controlling skewness.
#' @param tau Parameter controlling tailweight (kurtosis); must be positive (on log scale).
#'
#' @return A numeric vector of transformed values.
#' @keywords internal

#' @describeIn SHASH_transform Transforms SHASH-distributed data to approximately normal data.
#' @export
SHASH_to_normal <- function(x, mu, sigma, nu, tau){
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(mu) && length(mu) == 1)
  stopifnot(is.numeric(sigma) && length(sigma) == 1)
  stopifnot(is.numeric(nu) && length(nu) == 1)
  stopifnot(is.numeric(tau) && length(tau) == 1)
  
  sigma <- exp(sigma)
  tau   <- exp(tau)
  
  return(sinh((tau * asinh((x - mu)/ (sigma * tau))) - nu))
}

#' @describeIn SHASH_transform Transforms standard normal data back to SHASH-distributed scale.
#' @export
normal_to_SHASH <- function(x, mu, sigma, nu, tau){
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(mu) && length(mu) == 1)
  stopifnot(is.numeric(sigma) && length(sigma) == 1)
  stopifnot(is.numeric(nu) && length(nu) == 1)
  stopifnot(is.numeric(tau) && length(tau) == 1)
  
  sigma <- exp(sigma)
  tau   <- exp(tau)
  
  return((sigma * tau * sinh((asinh(x) + nu) / tau)) + mu)
}

#' SHASH-based Outlier Detection (Extended)
#'
#' Detect univariate outliers using an iterative SHASH fitting process with optional pre-flagging strategies.
#'
#' @param x Numeric vector.
#' @param thr0 Threshold for iterative convergence (default: 2.58).
#' @param thr Final outlier threshold (default: 4).
#' @param symmetric Logical. Use symmetric bounds for empirical rule.
#' @param use_huber Logical. Use Huber loss in empirical rule.
#' @param upper_only Logical. Flag only right-tail outliers.
#' @param use_isotree Logical. Use isolation forest scores.
#' @param use_isoplus Logical. Use isolation forest with directional split.
#' @param thr_isotree Score threshold for isolation forest.
#' @param maxit Max iterations.
#' @param weight_init Optional logical vector of initial weights.
#'
#' @return A list of class "SHASH_out".
#' @export
SHASH_out <- function(x,
                      thr0         = 2.58,
                      thr          = 4,
                      symmetric    = TRUE,
                      use_huber    = FALSE,
                      upper_only   = FALSE,
                      use_isotree  = FALSE,
                      use_isoplus  = FALSE,
                      thr_isotree  = 0.55,
                      maxit        = 100,
                      weight_init  = NULL) {
  if ((use_isotree || use_isoplus) && (thr_isotree < 0 || thr_isotree > 1)) {
    warning("`thr_isotree` should be between 0 and 1 when using isolation forest.")
  }
  if (use_isoplus && use_isotree) {
    stop("Cannot use `use_isoplus = TRUE` when `use_isotree = TRUE`. Set `use_isotree = FALSE`.")
  }
  if (!is.numeric(thr_isotree) || length(thr_isotree) != 1 || thr_isotree < 0 || thr_isotree > 1) stop("Threshold Isolation Forest 'thr_isotree' must be a positive numeric value between 0 and 1.")
  
  params <- list(
    orig_values  = x,
    thr0         = thr0,
    thr          = thr,
    symmetric    = symmetric,
    use_huber    = use_huber,
    upper_only   = upper_only,
    use_isotree  = use_isotree,
    use_isoplus  = use_isoplus,
    thr_isotree  = thr_isotree,
    maxit        = maxit,
    weight_init  = weight_init
  )
  
  na_locs   <- is.na(x)
  x_clean   <- x[!na_locs]
  n_clean   <- length(x_clean)
  x_med     <- median(x_clean)
  
  if (!is.null(weight_init)) {
    weight_new <- as.logical(weight_init[!na_locs])
    iso_scores <- rep(NA_real_, n_clean)
  } else if (use_isoplus) {
    iso_model  <- isolation.forest(data.frame(x_clean), ntrees=200)
    iso_scores <- predict(iso_model, data.frame(x_clean), type="score")
    idx_out0   <- which(iso_scores >= thr_isotree)
    if (length(idx_out0) == 0) {
      weight_new <- rep(TRUE, n_clean)
    } else {
      upp <- idx_out0[x_clean[idx_out0] >= x_med]
      low <- idx_out0[x_clean[idx_out0] <= x_med]
      if (upper_only) {
        weight_new <- if (length(upp) == 0) rep(TRUE, n_clean) else x_clean < min(x_clean[upp], na.rm=TRUE)
      } else {
        cutoff_upper <- if (length(upp) == 0) Inf else min(x_clean[upp], na.rm=TRUE)
        cutoff_lower <- if (length(low) == 0) -Inf else max(x_clean[low], na.rm=TRUE)
        weight_new <- (x_clean > cutoff_lower) & (x_clean < cutoff_upper)
      }
    }
  } else if (use_isotree) {
    iso_model  <- isolation.forest(data.frame(x_clean), ntrees=200)
    iso_scores <- predict(iso_model, data.frame(x_clean), type="score")
    weight_new <- iso_scores <= thr_isotree
  } else {
    # --- Default Initialization Using Robust Empirical Rule ---
    # Inverts output of emprule_rob() so TRUE = inlier
    iso_scores <- rep(NA_real_, n_clean)
    W <- tryCatch(
      1 - emprule_rob(x_clean, thr = thr0,
                      symmetric = symmetric,
                      use_huber = use_huber,
                      upper_only = upper_only),
      error = function(e) rep(TRUE, n_clean)
    )
    weight_new <- as.logical(W)
  }

  initial_weights_clean <- weight_new
  
  norm_iters_clean <- matrix(NA_real_,    nrow=n_clean, ncol=maxit)
  indx_iters_clean <- matrix(NA_integer_, nrow=n_clean, ncol=maxit)
  iter <- 0
  success <- FALSE
  
  repeat {
    iter <- iter + 1
    weight_old <- weight_new
    mod <- gamlss::gamlssML(
      x_clean ~ 1, family="SHASHo2", maxit=1e4,
      weights=as.numeric(weight_new)
    )
    est <- gamlss::coefAll(mod)
    
    x_norm_clean <- SHASH_to_normal(
      x_clean, mu=est$mu, sigma=est$sigma,
      nu=est$nu, tau=est$tau
    )
    norm_iters_clean[, iter] <- x_norm_clean
    
    weight_new <- if (upper_only) {
      x_norm_clean <= thr0
    } else {
      abs(x_norm_clean) <= thr0
    }
    indx_iters_clean[, iter] <- as.integer(!weight_new)
    
    if (isTRUE(all.equal(weight_old, weight_new))) {
      success <- TRUE
      break
    }
    if (iter >= maxit) break
  }
  
  norm_iters_clean <- norm_iters_clean[, seq_len(iter), drop=FALSE]
  indx_iters_clean <- indx_iters_clean[, seq_len(iter), drop=FALSE]
  
  final_idx_clean <- if (upper_only) {
    which(x_norm_clean > thr)
  } else {
    which(abs(x_norm_clean) > thr)
  }
  
  n_orig <- length(x)
  pad_vec <- function(v) {
    v_full <- rep(NA, n_orig)
    v_full[!na_locs] <- v
    v_full
  }
  x_norm_full         <- pad_vec(x_norm_clean)
  initial_weights_full<- pad_vec(initial_weights_clean)
  isotree_scores_full <- pad_vec(iso_scores)
  norm_iters_full     <- matrix(NA_real_,    nrow=n_orig, ncol=iter)
  norm_iters_full[!na_locs, ] <- norm_iters_clean
  indx_iters_full     <- matrix(NA_integer_, nrow=n_orig, ncol=iter)
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
  return(out)
}

#' Robust Empirical Rule Outlier Detection
#'
#' Detects outliers using the median ± threshold × MAD rule.
#'
#' @param x Numeric vector.
#' @param thr Threshold multiplier (default = 4).
#' @param symmetric Logical. If FALSE, use only upper tail.
#' @param use_huber Logical. If TRUE, apply Huber-like soft rejection.
#' @param upper_only Logical. If TRUE, use right-tail only.
#'
#' @return Logical vector: TRUE = outlier, FALSE = inlier.
#' @keywords internal
emprule_rob <- function(x, thr = 4, symmetric = TRUE, use_huber = FALSE, upper_only = FALSE) {
  x_med <- median(x, na.rm = TRUE)
  mad_val <- 1.4826 * median(abs(x - x_med), na.rm = TRUE)
  
  z <- (x - x_med) / mad_val
  if (use_huber) {
    z <- pmin(pmax(z, -thr), thr)  # shrink extreme z-scores
  }
  
  if (upper_only) {
    out <- z > thr
  } else if (symmetric) {
    out <- abs(z) > thr
  } else {
    out <- z < -thr | z > thr
  }
  
  return(out)
}


