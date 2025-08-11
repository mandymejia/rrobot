#' Plot Method for RD Analysis Results
#'
#' Creates diagnostic plots for robust distance analysis results.
#'
#' @param x An object of class "RD" from RD() or threshold_RD().
#' @param method Character string specifying threshold method. Auto-detected if NULL.
#' @param ... Additional arguments passed to method-specific plotting functions.
#' @return A ggplot object.
#' @method plot RD
#' @export
plot.RD <- function(x, method = NULL, ...) {
  # Auto-detect method if not specified
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

  # Dispatch to method-specific plotting
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
#' @param RD_obj RD_result object from compute_RD().
#' @param alpha Significance level for threshold label (default = 0.01).
#' @param binwidth Histogram bin width (default = 0.1).
#' @param log Logical. Apply log10 transformation to x-axis (default = TRUE).
#' @param show_f_density Logical. Show F-distribution curve overlay (default = TRUE).
#' @param ... Additional arguments.
#' @return A ggplot object.
#' @export
plot_F_histogram <- function(F_result, RD_obj, alpha = 0.01, binwidth = 0.1, log = FALSE, show_f_density = TRUE, ...) {
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
    (if (log) scale_x_log10() else scale_x_continuous(n.breaks = 10)) +
    ylab("Density") +
    xlab(if (log) "log10(Squared robust distances)" else "Squared robust distances") +
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
