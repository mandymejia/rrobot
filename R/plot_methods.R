#' Plot Method for RD Analysis Results
#'
#' Creates diagnostic plots for robust distance analysis results.
#'
#' @param x An object of class "RD" from RD() or threshold_RD().
#' @param type Character string specifying plot type: "histogram" (default), "imputations", or "univOut".
#' @param method Character string specifying threshold method. Auto-detected if NULL.
#' @param ... Additional arguments passed to plotting functions.
#' @return A ggplot object for \code{type = "histogram"} and \code{type = "univOut"};
#'   for \code{type = "imputations"} the plots are printed and returned invisibly as a list.
#' @method plot RD
#' @export
plot.RD <- function(x, type = c("histogram", "imputations", "univOut"), method = NULL, ...) {
  type <- match.arg(type)
  if (type == "univOut") {
    return(plot_univOut(x, ...))
  }

  if (type == "imputations") {
    return(plot_imputations(x, ...))
  }

  # type == "histogram" (default)
  th <- x$thresholds
  is_multi <- identical(class(th), "list")  # method = "all": plain list of results

  if (is.null(method)) {
    if (is_multi) {
      return(plot_RD_histogram_multi(x, x$RD_obj, ...))
    }
    method <- detect_threshold_method(th)
  } else {
    method <- match.arg(method, c("SI", "SI_boot", "MI", "MI_boot", "F", "SHASH"))
  }

  thresh_obj <- if (is_multi) th[[method]] else th
  if (is.null(thresh_obj)) {
    stop("Threshold method '", method, "' is not available in this RD object.")
  }

  # SHASH thresholds live on the distances of the SHASH-normalised data
  rd_obj <- if (method == "SHASH" && !is.null(x$RD_obj_shash)) x$RD_obj_shash else x$RD_obj

  if (method == "F") {
    plot_F_histogram(thresh_obj, rd_obj, ...)
  } else {
    plot_RD_histogram(thresh_obj, rd_obj, ...)
  }
}

# Map the class of a threshold result to its method name.
#' @noRd
detect_threshold_method <- function(th) {
  map <- c(F_result = "F", SI_result = "SI", SI_boot_result = "SI_boot",
           MI_result = "MI", MI_boot_result = "MI_boot", SHASH_result = "SHASH")
  hit <- intersect(names(map), class(th))
  if (length(hit) == 0) stop("Unknown threshold result type.")
  unname(map[hit[1]])
}

# Included/excluded status of each observation as a factor.
#' @noRd
rd_obs_status <- function(n, ind_incld) {
  factor(ifelse(seq_len(n) %in% ind_incld, "included", "excluded"),
         levels = c("included", "excluded"))
}

# Shared histogram + threshold line used by plot_RD_histogram() and plot_F_histogram().
#' @noRd
RD_hist_base <- function(df_RD, q_HR, binwidth) {
  ggplot(df_RD, aes(x = .data$distance)) +
    geom_histogram(aes(y = after_stat(density), fill = .data$observation, color = .data$observation),
                   binwidth = binwidth, alpha = 0.5) +
    geom_vline(xintercept = q_HR, linetype = "dashed", linewidth = 0.8, color = "#333333") +
    scale_fill_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
    scale_color_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
    scale_x_log10() +
    ylab("Density") +
    xlab("log10(Squared robust distances)") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 13),
      strip.text.y = element_text(angle = 0),
      legend.key = element_blank()
    )
}

# Use excluded distances when possible; fall back to all distances for sparse
# exclusions so density estimation remains defined.
#' @noRd
plot_peak_density <- function(values, fallback_values) {
  values <- values[is.finite(values)]
  if (length(values) < 2L) {
    values <- fallback_values[is.finite(fallback_values)]
  }
  if (length(values) < 2L) return(0)

  max(stats::density(values, from = 0)$y)
}

#' Plot Multiple Threshold Methods on Robust Distance Histogram
#'
#' Creates a histogram of robust distances with multiple colored threshold lines
#' showing different outlier detection methods simultaneously.
#'
#' @param RD_result An RD result object from threshold_RD() with method="all" containing
#'   multiple threshold results in a list.
#' @inheritParams params_doc
#' @param methods Character vector of threshold methods to display (default: c("SI", "SI_boot", "MI", "MI_boot", "F")).
#'
#' @return A ggplot object with histogram colored by inclusion status and multiple
#'   colored threshold lines for comparison of different methods.
#' @keywords internal
plot_RD_histogram_multi <- function(RD_result, RD_obj, methods = c("SI", "SI_boot", "MI", "MI_boot", "F"), ...) {
  stopifnot("RD_result$thresholds must be a list of threshold results (use method = 'all')" =
              identical(class(RD_result$thresholds), "list"))
  stopifnot("RD_obj must have $RD and $ind_incld elements" = !is.null(RD_obj$RD) && !is.null(RD_obj$ind_incld))

  method_labels <- c(
    "SI" = "SI (99%)",
    "SI_boot" = "SI_boot (2.5% CI LB)",
    "MI" = "MI (2.5%)",
    "MI_boot" = "MI_boot",
    "F" = "F"
  )
  stopifnot("methods must be a subset of SI, SI_boot, MI, MI_boot, F" = all(methods %in% names(method_labels)))
  stopifnot("Requested methods are not available in RD_result$thresholds" =
              all(methods %in% names(RD_result$thresholds)))

  # --- Data for histogram -------------------------------------------------------
  RD      <- RD_obj$RD
  eps     <- 1e-12
  df_RD   <- data.frame(z = log(pmax(RD, eps)),
                        label = rd_obs_status(length(RD), RD_obj$ind_incld))

  # --- Collect thresholds (on RD scale) -----------------------------------------
  # Every method's $threshold is already on the RD scale (for "F" it is the
  # F quantile divided by the scale factor), so no rescaling is needed.
  thresholds <- vapply(methods, function(m) RD_result$thresholds[[m]]$threshold, numeric(1))

  df_thresh <- data.frame(
    method  = factor(method_labels[methods], levels = method_labels[methods]),
    thr_log = log(pmax(thresholds, eps))
  )

  # --- Plot ---------------------------------------------------------------------
  ggplot(df_RD, aes(x = .data$z, fill = .data$label)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 30, position = "identity", alpha = 0.6, color = "black") +
    geom_vline(data = df_thresh,
               aes(xintercept = .data$thr_log, color = .data$method),
               linewidth = 0.9, linetype = "solid", inherit.aes = FALSE) +
    scale_fill_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
    labs(
      title = "Histogram of log(RD) with Thresholds (non-SHASH)",
      x = "log(RD)",
      y = "Density",
      fill = "Observation",
      color = "Threshold"
    ) +
    theme_minimal(base_size = 12)
}

#' Plot Robust Distance Histogram with Threshold
#'
#' Creates a histogram of robust distances with threshold line for outlier detection.
#'
#' @inheritParams params_doc
#'
#' @return A ggplot object with histogram colored by inclusion status and threshold line.
#' @keywords internal
plot_RD_histogram <- function(thresh_result, RD_obj, binwidth = 0.1, ...) {

  stopifnot("thresh_result must have a $threshold element" = !is.null(thresh_result$threshold))
  stopifnot("RD_obj must have $RD and $ind_incld elements" = !is.null(RD_obj$RD) && !is.null(RD_obj$ind_incld))

  RD   <- RD_obj$RD
  q_HR <- thresh_result$threshold

  df_RD <- data.frame(
    distance = RD,
    observation = rd_obs_status(length(RD), RD_obj$ind_incld)
  )
  df_RD <- df_RD[df_RD$distance > 0, ]  # log10 axis

  RD_hist_base(df_RD, q_HR, binwidth)
}


#' Plot F-Distribution Method Results
#'
#' Creates histogram of robust distances with F-distribution overlay and threshold.
#'
#' @param F_result F_result object from thresh_F().
#' @inheritParams params_doc
#' @param show_f_density Logical. Show F-distribution curve overlay (default = TRUE).
#' @return A ggplot object.
#' @keywords internal
plot_F_histogram <- function(F_result, RD_obj, binwidth = 0.1, show_f_density = TRUE, ...) {
  stopifnot("RD_obj must have $RD and $ind_incld elements" = !is.null(RD_obj$RD) && !is.null(RD_obj$ind_incld))

  RD <- RD_obj$RD

  # Get threshold and scaling from F_result
  scale_factor <- F_result$scale
  q_HR <- F_result$threshold * scale_factor
  df1 <- F_result$df[1]
  df2 <- F_result$df[2]

  observation <- rd_obs_status(length(RD), RD_obj$ind_incld)
  df_RD <- data.frame(distance = RD * scale_factor, observation = observation)

  p <- RD_hist_base(df_RD[df_RD$distance > 0, ], q_HR, binwidth)

  if (show_f_density) {
    # Extend tail range
    excl_vals <- RD[observation == "excluded"] * scale_factor
    peak_density <- plot_peak_density(excl_vals, RD * scale_factor)
    xvals <- seq(0, max(c(excl_vals, q_HR), na.rm = TRUE) * 1.1, length.out = 512)
    xvals <- xvals[xvals > 0]  # log10 axis

    # F-distribution curve scaled to match empirical peak
    f_density <- stats::df(xvals, df1 = df1, df2 = df2)
    yvals <- f_density * (peak_density / max(f_density))

    p <- p + geom_line(data = data.frame(x = xvals, y = yvals),
                       aes(x = .data$x, y = .data$y),
                       linewidth = 0.8, color = "#333333", inherit.aes = FALSE)
  }

  p
}

#' Plot Univariate Outliers from RD Analysis
#'
#' Creates a heatmap visualization of univariate outliers detected in high-kurtosis components.
#'
#' @param x An object of class "RD" from RD() or threshold_RD().
#' @inheritParams params_doc
#'
#' @return A ggplot object showing a heatmap of outlier locations.
#' @keywords internal
plot_univOut <- function(x, ...) {
  imp_x <- attr(x, "univOut_hk")
  stopifnot("univOut data not available" = !is.null(imp_x))

  out_hk <- imp_x$outliers

  # ---- Heatmap of hk outlier matrix (long format, no reshape2 needed) ----
  mat_df <- data.frame(
    Time     = rep(seq_len(nrow(out_hk)), times = ncol(out_hk)),
    Variable = rep(seq_len(ncol(out_hk)), each = nrow(out_hk)),
    Outlier  = as.vector(out_hk)
  )

  ggplot2::ggplot(mat_df, ggplot2::aes(x = .data$Variable, y = .data$Time, fill = .data$Outlier)) +
    ggplot2::geom_tile(color = "#D6E9FF") +
    ggplot2::scale_fill_manual(values = c("FALSE" = "#f5f5f5", "TRUE" = "#191970")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Outlier Matrix", x = "Variable", y = "Observations") +
    ggplot2::theme(aspect.ratio = 1.2)
}

#' Plot Multiple Imputation Results from RD Analysis
#'
#' Creates time series plots showing original data, temporal imputation,
#' and multiple imputation results with outlier locations highlighted.
#'
#' @param x An object of class "RD" from RD() or threshold_RD().
#' @inheritParams params_doc
#'
#' @return Prints one ggplot object per variable and returns them invisibly as a list.
#' @keywords internal
plot_imputations <- function(x, ...) {
  x_data       <- attr(x, "x_data")
  out_result   <- attr(x, "univOut_hk")
  imp_result   <- attr(x, "impute_univOut_hk")
  multiple_imp <- attr(x, "MImpute_hk")
  stopifnot("Multiple imputation data not available. This plot only available for MI and MI_boot methods." =
              !is.null(multiple_imp))

  x_data <- as.matrix(x_data)
  n_time <- nrow(x_data)
  NA_locs <- imp_result$NA_locs
  imputed_datasets <- multiple_imp$imp_datasets
  num_imputations <- length(imputed_datasets)
  num_columns <- ncol(out_result$outliers)  # Number of time series (columns)

  # Choose a random imputed dataset to highlight in orange
  random_imputed_idx <- sample.int(num_imputations, 1)

  # Imputation method used (stored by impute_univOut(); older objects fall back to the call)
  method_used <- imp_result$method
  if (is.null(method_used)) {
    method_used <- if ("interp" %in% as.character(x$call)) "interp" else "mean"
  }
  message("Using method: ", method_used)

  plots <- vector("list", num_columns)

  # Loop through each column (IC)
  for (q in seq_len(num_columns)) {
    plot_data <- data.frame(Time = seq_len(n_time), Original = x_data[, q])
    # NA (outlier) locations for this column
    na_rows <- NA_locs[NA_locs[, 2] == q, 1]

    # All imputed datasets for this column in one long data frame (one layer instead of M)
    mi_long <- data.frame(
      Time    = rep(seq_len(n_time), times = num_imputations),
      Imputed = unlist(lapply(imputed_datasets, function(d) d[, q]), use.names = FALSE),
      imp     = rep(seq_len(num_imputations), each = n_time)
    )
    is_random <- mi_long$imp == random_imputed_idx

    p <- ggplot() +
      # red dots for the outlier (now imputed) locations
      geom_point(data = data.frame(Time = na_rows, Original = x_data[na_rows, q]),
                 aes(x = .data$Time, y = .data$Original, color = "Outliers"), size = 2, alpha = 0.7) +
      # temporal imputation line
      geom_line(data = data.frame(Time = seq_len(n_time), Temporal = imp_result$imp_data[, q]),
                aes(x = .data$Time, y = .data$Temporal, color = "Temporal Imputation"),
                linewidth = 1.0, linetype = "dashed") +
      # imputed datasets as transparent lines
      geom_line(data = mi_long[!is_random, ],
                aes(x = .data$Time, y = .data$Imputed, group = .data$imp, color = "Multiple Imputations"),
                linewidth = 0.5, alpha = 0.02) +
      # the randomly selected imputation
      geom_line(data = mi_long[is_random, ],
                aes(x = .data$Time, y = .data$Imputed, color = "Random MI"),
                linewidth = 0.8)

    # Original time series (plotted last)
    if (method_used == "mean") {
      p <- p + geom_point(data = plot_data, aes(x = .data$Time, y = .data$Original, color = "Original"), size = 0.8)
    } else {
      p <- p + geom_line(data = plot_data, aes(x = .data$Time, y = .data$Original, color = "Original"), linewidth = 0.8)
    }

    p <- p + scale_color_manual(
      name = "Data Type",
      values = c(
        "Original" = "black",
        "Temporal Imputation" = "green",
        "Multiple Imputations" = "blue",
        "Random MI" = "orange",
        "Outliers" = "red"
      )
    ) +
      labs(title = paste("Time Series for Column", q, "with Missing Data and Imputed Values"),
           x = "Time", y = "Value") +
      theme_minimal() +
      theme(legend.position = "bottom")

    print(p)
    plots[[q]] <- p
  }

  invisible(plots)
}
