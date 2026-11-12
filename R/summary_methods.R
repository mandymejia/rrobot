#' Summary method for SI results
#' @param object An object of class "SI_result"
#' @inheritParams params_doc
#' @method summary SI_result
#' @return NULL, invisibly
#' @export
summary.SI_result <- function(object, ...) {
  cat("Single Imputation (SI) Results\n")
  cat("==============================\n")
  cat("SI Threshold:", round(object$threshold, 4), "\n")
  cat("Number of observations:", length(object$SI_obj$RD), "\n")
  cat("Robust subset size:", object$SI_obj$h, "\n")

  # Outliers are the original observations whose RD exceeds the threshold
  n_outliers <- length(object$flagged_outliers)
  pct_outliers <- round(100 * n_outliers / length(object$SI_obj$RD), 1)
  cat("Outliers detected:", paste0(n_outliers, " (", pct_outliers, "%)"), "\n")
  cat("\n")

  invisible(NULL)
}

#' Summary method for SI_boot results
#' @param object An object of class "SI_boot_result"
#' @inheritParams params_doc
#' @method summary SI_boot_result
#' @return NULL, invisibly
#' @export
summary.SI_boot_result <- function(object, ...) {
  cat("Single Imputation Bootstrap (SI_boot) Results\n")
  cat("=============================================\n")
  cat("Bootstrap samples:", length(object$thresholds), "\n")
  cat("Confidence Interval:\n")
  cat("  Lower bound (threshold):", round(object$threshold, 4), "\n")
  cat("  Upper bound:", round(object$UB_CI, 4), "\n")
  cat("Flagged outliers:", length(object$flagged_outliers), "\n")
  cat("\n")

  invisible(NULL)
}

#' Summary method for MI results
#' @param object An object of class "MI_result"
#' @inheritParams params_doc
#' @method summary MI_result
#' @return NULL, invisibly
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
  cat("Flagged outliers:", length(object$flagged_outliers), "\n")
  cat("\n")
  invisible(NULL)
}

#' Summary method for MI_boot results
#' @param object An object of class "MI_boot_result"
#' @inheritParams params_doc
#' @method summary MI_boot_result
#' @return NULL, invisibly
#' @export
summary.MI_boot_result <- function(object, ...) {
  cat("Multiple Imputation Bootstrap (MI_boot) Results\n")
  cat("===============================================\n")

  total_bootstraps <- length(object$thresholds)
  cat("Total bootstrap samples:", total_bootstraps, "\n")
  cat("Final threshold:", round(object$threshold, 4), "\n")

  # flagged_outliers holds row indices, so count them (do not sum them)
  n_outliers <- length(object$flagged_outliers)
  if (!is.null(object$n_obs)) {
    pct_outliers <- round(100 * n_outliers / object$n_obs, 1)
    cat("Flagged outliers:", paste0(n_outliers, " (", pct_outliers, "%)"), "\n")
  } else {
    cat("Flagged outliers:", n_outliers, "\n")
  }
  cat("\n")

  invisible(NULL)
}

# Shared body for the F and SHASH summaries
summarize_F_result <- function(object, title, rule = strrep("=", nchar(title))) {
  cat(title, "\n", sep = "")
  cat(rule, "\n", sep = "")
  cat("F-distribution threshold:", round(object$threshold, 4), "\n")
  cat("Degrees of freedom:", paste(round(object$df, 2), collapse = ", "), "\n")
  cat("Consistency factor (c):", round(object$c, 4), "\n")
  cat("Scale factor:", round(object$scale, 4), "\n")
  cat("Estimated m parameter:", round(object$m, 2), "\n")
  cat("Flagged outliers:", length(object$flagged_outliers), "\n")
  cat("\n")

  invisible(NULL)
}

#' Summary method for Hardin & Rocke F results
#' @param object An object of class "F_result"
#' @inheritParams params_doc
#' @method summary F_result
#' @return NULL, invisibly
#' @export
summary.F_result <- function(object, ...) {
  summarize_F_result(object, "Hardin & Rocke F-Distribution Results", rule = "====================================")
}

#' Summary method for SHASH results
#' @param object An object of class "SHASH_result"
#' @inheritParams params_doc
#' @method summary SHASH_result
#' @return NULL, invisibly
#' @export
summary.SHASH_result <- function(object, ...) {
  summarize_F_result(object, "SHASH (F-Distribution Threshold) Results")
}

#' Print method for RD results
#'
#' Compact display of an \code{"RD"} object (the stored imputations and data
#' attached to it are not printed).
#'
#' @param x An object of class "RD" from threshold_RD().
#' @inheritParams params_doc
#' @method print RD
#' @return \code{x}, invisibly
#' @export
print.RD <- function(x, ...) {
  cat("Robust distance outlier detection (class 'RD')\n")
  if (!is.null(x$RD_obj$RD)) {
    cat("Observations:", length(x$RD_obj$RD), "\n")
  }

  describe <- function(th) {
    paste0("threshold = ", format(signif(th$threshold, 5)),
           ", flagged = ", length(th$flagged_outliers))
  }

  th <- x$thresholds
  if (identical(class(th), "list")) {
    # method = "all": a named list of threshold results
    for (nm in names(th)) cat(sprintf("  %-8s %s\n", nm, describe(th[[nm]])))
  } else {
    cat(sprintf("  %-8s %s\n", sub("_result$", "", class(th)[1]), describe(th)))
  }

  invisible(x)
}
