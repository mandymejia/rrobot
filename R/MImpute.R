#' Multiple Imputation with Per-Cycle Updates (OLS + MICE-style)
#'
#' @param x  (T × p) high-kurtosis ICA matrix to impute.
#' @param w  (T × q) predictors (e.g., low-kurtosis components).
#' @param outlier_matrix logical (T × p) mask of entries to impute.
#' @param M  number of multiply-imputed datasets (default 50).
#' @param k  number of chained-equation cycles per dataset (default 5–10 is common).
#' @param ridge_eps tiny ridge added to X'X for stability (default 1e-8).
#' @param tol optional early-stop tolerance on per-cycle max change (NA to disable).
#' @return list(imp_datasets, outlier_matrix)
#' @importFrom MASS mvrnorm ginv
#' @importFrom stats var rnorm
#' @export
MImpute <- function(x, w, outlier_matrix, M = 50, k = 5, ridge_eps = 1e-8, tol = NA_real_) {
  stopifnot(is.matrix(x), is.matrix(w), is.logical(outlier_matrix))
  stopifnot(all(dim(x) == dim(outlier_matrix)))
  stopifnot(nrow(w) == nrow(x))

  Tn <- nrow(x)
  p  <- ncol(x)
  imp_datasets <- vector("list", M)

  for (m in seq_len(M)) {
    x_imp <- x

    for (cyc in seq_len(k)) {
      max_change <- 0

      for (j in seq_len(p)) {
        idx_imp <- which(outlier_matrix[, j])
        if (length(idx_imp) == 0) next

        yj <- x_imp[, j]
        yj[idx_imp] <- NA

        obs <- !is.na(yj)
        if (sum(obs) < 10) next

        # Predictors use current imputations for other hk-data vars
        if (p > 1) {
          other_x <- x_imp[, -j, drop = FALSE]
          P <- cbind(w, other_x)
        } else {
          P <- w
        }

        # Design matrices
        Xobs <- cbind(Intercept = 1, P[obs, , drop = FALSE])
        yobs <- yj[obs]
        Xnew <- cbind(Intercept = 1, P[idx_imp, , drop = FALSE])

        # OLS fit with tiny ridge
        XtX <- crossprod(Xobs)
        Xty <- crossprod(Xobs, yobs)
        XtX_ridge <- XtX + diag(ridge_eps, ncol(XtX))
        beta_hat <- tryCatch(
          solve(XtX_ridge, Xty),
          error = function(e) MASS::ginv(XtX_ridge) %*% Xty
        )

        # Residual variance and beta covariance
        qrX <- qr(Xobs); rnk <- qrX$rank
        resid  <- as.numeric(yobs - Xobs %*% beta_hat)
        sigma2 <- if (length(yobs) > rnk) sum(resid^2) / (length(yobs) - rnk) else var(resid)
        XtX_inv <- tryCatch(solve(XtX), error = function(e) MASS::ginv(XtX))
        beta_cov <- sigma2 * (XtX_inv + diag(ridge_eps, ncol(XtX_inv)))  # PD for mvrnorm

        # Proper MI: one beta draw + residual noise
        b_draw <- tryCatch(
          MASS::mvrnorm(1, mu = as.numeric(beta_hat), Sigma = beta_cov),
          error = function(e) as.numeric(beta_hat)
        )
        mu_hat <- as.numeric(Xnew %*% b_draw)
        y_draw <- mu_hat + rnorm(length(idx_imp), sd = sqrt(max(sigma2, 1e-12)))

        # Update immediately
        old_vals <- x_imp[idx_imp, j]
        x_imp[idx_imp, j] <- y_draw
        max_change <- max(max_change, max(abs(y_draw - old_vals), na.rm = TRUE))
      }
    }

    imp_datasets[[m]] <- x_imp
  }

  list(imp_datasets = imp_datasets,
       outlier_matrix = outlier_matrix)
}
