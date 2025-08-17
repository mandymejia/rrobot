#' Plot Method for RD Analysis Results
#'
#' Creates diagnostic plots for robust distance analysis results.
#'
#' @param x An object of class "RD" from RD() or threshold_RD().
#' @param type Character string specifying plot type: "histogram" (default), "thresholds", "imputation", or "univOut".
#' @param method Character string specifying threshold method. Auto-detected if NULL.
#' @param ... Additional arguments passed to plotting functions.
#' @return A ggplot object.
#' @method plot RD
#' @export
plot.RD <- function(x, type = c("histogram", "thresholds", "imputation", "univOut"), method = NULL, ...) {
  type <- match.arg(type)
  if (type == "univOut") {
    return(plot_univOut(x, ...))
  }

  if (type == "imputation") {
    return(plot_imputation(x, ...))
  }

  if (type == "thresholds") {
    return(plot_thresholds(x, ...))
  }

  # type == "histogram" (default) - auto-detect method and show histogram
  if (is.null(method)) {
    if (inherits(x$thresholds, "F_result")) {
      method <- "F"
    } else if (inherits(x$thresholds, "SI_result")) {
      method <- "SI"
    } else if (inherits(x$thresholds, "list")) {
      stop("Multiple methods detected. Please specify method parameter.")
    }
    # Add other method detection...
  }

  # Dispatch to method-specific plotting for default histograms
  switch(method,
         "F" = plot_F_histogram(x$thresholds, x$RD_obj, ...),
         "SI" = plot_SI_method(x$thresholds, x$RD_obj, ...),
         # Add other methods...
         stop("Unsupported method: ", method)
  )
}

#' Plot F-Distribution Method Results
#'
#' Creates histogram of robust distances with F-distribution overlay and threshold.
#'
#' @param F_result F_result object from thresh_F().
#' @inheritParams RD_obj
#' @inheritParams alpha
#' @param binwidth Histogram bin width (default = 0.1).
#' @param log Logical. Apply log10 transformation to x-axis (default = TRUE).
#' @param show_f_density Logical. Show F-distribution curve overlay (default = TRUE).
#' @inheritParams ...
#' @return A ggplot object.
#' @keywords internal
plot_F_histogram <- function(F_result, RD_obj, alpha = 0.01, binwidth = 0.1, show_f_density = TRUE, ...) {
  # Extract what we need from the pre-computed objects
  RD <- RD_obj$RD
  ind_incld <- RD_obj$ind_incld
  t <- length(RD)

  # Get threshold and scaling from F_result
  q_HR <- F_result$threshold * F_result$scale
  scale_factor <- F_result$scale
  df1 <- F_result$df[1]
  df2 <- F_result$df[2]

  # Create observation status
  observation <- ifelse(seq_len(t) %in% ind_incld, "included", "excluded")

  # Extend tail range
  excl_vals <- RD[observation == "excluded"] * scale_factor
  dens_excld <- density(excl_vals, from = 0)
  peak_density <- max(dens_excld$y)
  xvals <- seq(0, max(c(excl_vals, q_HR), na.rm = TRUE) * 1.1, length.out = 512)

  # F-distribution curve scaled to match empirical peak
  f_density <- df(xvals, df1 = df1, df2 = df2)
  scale_factor_y <- peak_density / max(f_density)
  yvals <- f_density * scale_factor_y

  # Create data frames for plotting (single dataset, no "data" column needed)
  df_RD <- data.frame(
    distance = RD * scale_factor,
    observation = observation
  )

  df_HR <- data.frame(x = xvals, y = yvals, q_HR = q_HR)
  df_HR_label <- data.frame(
    x = q_HR,
    y = max(yvals, na.rm = TRUE) * 1.05,
    label = paste0("F(", 1 - alpha, ")")
  )

  # Filter positive values
  df_RD <- df_RD[df_RD$distance > 0, ]
  df_HR <- df_HR[df_HR$q_HR > 0, ]
  df_HR_label <- df_HR_label[df_HR_label$x > 0, ]

  # Factor the observation column
  df_RD$observation <- factor(df_RD$observation, levels = c("included", "excluded"))

  # Create the plot (no faceting for single dataset)
  p <- ggplot(df_RD, aes(x = .data$distance)) +
    geom_histogram(aes(y = after_stat(density), fill = observation, color = observation),
                   binwidth = binwidth, alpha = 0.5) +
    geom_vline(xintercept = q_HR, linetype = "dashed", linewidth = 0.8, color = '#333333') +
    (if (show_f_density) geom_line(data = df_HR, aes(x = .data$x, y = .data$y), size = 0.8, color = '#333333', inherit.aes = FALSE) else NULL) +
    geom_text(data = df_HR_label,
              aes(x = .data$x, y = .data$y, label = .data$label),
              inherit.aes = FALSE,
              hjust = -0.1, vjust = -0.5, size = 4.5, color = "#333333") +
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

  return(p)
}

#' Plot Univariate Outliers from RD Analysis
#'
#' Creates a heatmap visualization of univariate outliers detected in high-kurtosis components.
#'
#' @param x An object of class "RD" from RD() or threshold_RD().
#' @inheritParams cutoff
#' @inheritParams method_univOut
#'
#' @return A ggplot object showing a heatmap of outlier locations.
#' @keywords internal
plot_univOut <- function(
    x,
    cutoff = NULL,
    method = NULL
) {
  method <- match.arg(method)

  # Get attributes from S3
  imp_x <- attr(x, "univOut_hk")
  stopifnot("univOut data not available" = !is.null(imp_x))

  # Get cutoff/method if user specified
  call_obj <- x$call
  cutoff <- if (is.null(cutoff)) (if ("cutoff" %in% names(call_obj)) eval(call_obj$cutoff) else 4) else cutoff
  method <- if (is.null(method)) (if ("trans" %in% names(call_obj)) as.character(eval(call_obj$trans)) else "SHASH") else method

  # Use the extracted parameters (or function defaults if user didn't specify in call)
  final_cutoff <- cutoff  # Function parameter takes precedence
  final_method <- method  # Function parameter takes precedence

  message("Plotting univOut with method=", method, ", cutoff=", cutoff)

  out_hk <- imp_x$outliers

  # True indices and timepoints
  idx_hk <- which(out_hk, arr.ind = TRUE)
  tp_hk  <- unique(idx_hk[, 1])

  # ---- Heatmap of hk outlier matrix ----
  mat_df <- reshape2::melt(out_hk)
  colnames(mat_df) <- c("Time", "Variable", "Outlier")
  p <- ggplot2::ggplot(mat_df, ggplot2::aes(x = Variable, y = Time, fill = Outlier)) +
    ggplot2::geom_tile(color = "#D6E9FF")+
    ggplot2::scale_fill_manual(values = c("FALSE" = "#f5f5f5", "TRUE" = "#191970")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Outlier Matrix (High-Kurtosis Components)", x = "Variable", y = "Time") +
    ggplot2::theme(aspect.ratio = 1.2) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))

  return(p)
}
