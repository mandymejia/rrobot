#' Summary method for SI results
#' @param object An object of class "SI_result"
#' @inheritParams ...
#' @method summary SI_result
#' @export
summary.SI_result <- function(object, ...) {
  cat("Single Imputation (SI) Results\n")
  cat("==============================\n")
  cat("SI Threshold:", round(object$threshold, 4), "\n")
  cat("Number of observations:", length(object$SI_obj$RD), "\n")
  cat("Robust subset size:", object$SI_obj$h, "\n")

  # Can only show outlier count if we have the threshold and RD
  outliers <- object$SI_obj$RD > object$threshold
  n_outliers <- sum(outliers, na.rm = TRUE)
  pct_outliers <- round(100 * n_outliers / length(outliers), 1)
  cat("Outliers detected:", paste0(n_outliers, " (", pct_outliers, "%)"), "\n")
  cat("\n")

  invisible(object)
}

#' Summary method for SI_boot results
#' @param object An object of class "SI_boot_result"
#' @inheritParams ...
#' @method summary SI_boot_result
#' @export
summary.SI_boot_result <- function(object, ...) {
  cat("Single Imputation Bootstrap (SI_boot) Results\n")
  cat("=============================================\n")
  cat("Bootstrap samples:", length(object$thresholds), "\n")
  cat("Confidence Interval:\n")
  cat("  Lower bound (threshold):", round(object$threshold, 4), "\n")
  cat("  Upper bound:", round(object$UB_CI, 4), "\n")
  cat("\n")

  invisible(object)
}

#' Summary method for MI results
#' @param object An object of class "MI_result"
#' @inheritParams ...
#' @method summary MI_result
#' @export
summary.MI_result <- function(object, ...) {
  cat("Multiple Imputation (MI) Results\n")
  cat("================================\n")
  cat("Number of imputations:", length(object$thresholds), "\n")
  cat("Threshold summary:\n")
  cat("  Min:", round(min(object$thresholds), 4), "\n")
  cat("  Mean:", round(mean(object$thresholds), 4), "\n")
  cat("  Max:", round(max(object$thresholds), 4), "\n")
  cat("Final threshold:", round(object$threshold, 4), "\n")
  cat("\n")
  invisible(object)
}

#' Summary method for MI_boot results
#' @param object An object of class "MI_boot_result"
#' @inheritParams ...
#' @method summary MI_boot_result
#' @export
summary.MI_boot_result <- function(object, ...) {
  cat("Multiple Imputation Bootstrap (MI_boot) Results\n")
  cat("===============================================\n")

  total_bootstraps <- length(object$thresholds)
  cat("Total bootstrap samples:", total_bootstraps, "\n")
  cat("Final threshold:", round(object$threshold, 4), "\n")

  n_outliers <- sum(object$flagged_outliers, na.rm = TRUE)
  pct_outliers <- round(100 * n_outliers / length(object$flagged_outliers), 1)

  cat("Flagged outliers:", paste0(n_outliers, " (", pct_outliers, "%)"), "\n")
  cat("\n")

  invisible(object)
}

#' Summary method for Hardin & Rocke F results
#' @param object An object of class "F_result" or "HR_result"
#' @inheritParams ...
#' @method summary F_result
#' @export
summary.F_result <- function(object, ...) {
  cat("Hardin & Rocke F-Distribution Results\n")
  cat("====================================\n")
  cat("F-distribution threshold:", round(object$threshold, 4), "\n")
  cat("Degrees of freedom:", paste(round(object$df, 2), collapse = ", "), "\n")
  cat("Consistency factor (c):", round(object$c, 4), "\n")
  cat("Scale factor:", round(object$scale, 4), "\n")
  cat("Estimated m parameter:", round(object$m, 2), "\n")
  cat("\n")

  invisible(object)
}

