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
plot.RD <- function(x, type = c("histogram", "thresholds", "imputations", "univOut"), method = NULL, ...) {
  type <- match.arg(type)
  if (type == "univOut") {
    return(plot_univOut(x, ...))
  }

  if (type == "imputations") {
    return(plot_imputations(x, ...))
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
    } else if (inherits(x$thresholds, "SI_boot_result")) {
      method <- "SI_boot"
    } else if (inherits(x$thresholds, "MI_result")) {
      method <- "MI"
    } else if (inherits(x$thresholds, "MI_boot_result")) {
      method <- "MI_boot"
    } else if (inherits(x$thresholds, "SHASH_result")) {
      method <- "SHASH"
    } else if (inherits(x$thresholds, "list")) {
      return(plot_RD_histogram_multi(x, x$RD_obj, ...))
    } else {
      stop("Unknown threshold result type.")
    }
  }

  thresh_obj <- if (inherits(x$thresholds, "list")) x$thresholds[[method]] else x$thresholds

  # Dispatch to method-specific plotting for default histograms
  switch(method,
         "F" = plot_F_histogram(thresh_obj, x$RD_obj, ...),
         "SI" = plot_RD_histogram(thresh_obj, x$RD_obj, ...),
         "SI_boot" = plot_RD_histogram(thresh_obj, x$RD_obj, ...),
         "MI" = plot_RD_histogram(thresh_obj, x$RD_obj, ...),
         "MI_boot" = plot_RD_histogram(thresh_obj, x$RD_obj, ...),
         "SHASH" = plot_RD_histogram(thresh_obj, x$RD_obj, ...),
         stop("Unsupported method: ", method)
  )
}

#' Plot Multiple Threshold Methods on Robust Distance Histogram
#'
#' Creates a histogram of robust distances with multiple colored threshold lines
#' showing different outlier detection methods simultaneously.
#'
#' @param RD_result An RD result object from threshold_RD() with method="all" containing
#'   multiple threshold results in a list.
#' @inheritParams RD_obj
#' @param methods Character vector of threshold methods to display (default: c("SI", "SI_boot", "MI", "MI_boot")).
#' @inheritParams alpha
#' @inheritParams binwidth
#' @inheritParams ...
#'
#' @return A ggplot object with histogram colored by inclusion status and multiple
#'   colored threshold lines for comparison of different methods.
#' @keywords internal
plot_RD_histogram_multi <- function(RD_result, RD_obj, methods = c("SI", "SI_boot", "MI", "MI_boot"), alpha = 0.01, binwidth = 0.1, ...) {
  stopifnot("thresh_result must have a $threshold element" = inherits(RD_result$thresholds, "list"))
  stopifnot("RD_obj must have $RD and $ind_incld elements" = !is.null(RD_obj$RD) && !is.null(RD_obj$ind_incld))

  # Extract what we need from the pre-computed objects
  RD <- RD_obj$RD
  ind_incld <- RD_obj$ind_incld
  t <- length(RD)

  scale_factor <- 1

  # Create observation status
  observation <- ifelse(seq_len(t) %in% ind_incld, "included", "excluded")

  # Extend tail range
  excl_vals <- RD[observation == "excluded"] * scale_factor
  dens_excld <- density(excl_vals, from = 0)
  peak_density <- max(dens_excld$y)

  # Create data frames for plotting (single dataset, no "data" column needed)
  df_RD <- data.frame(
    distance = RD * scale_factor,
    observation = observation
  )

  # Filter positive values
  df_RD <- df_RD[df_RD$distance > 0, ]

  # Factor the observation column
  df_RD$observation <- factor(df_RD$observation, levels = c("included", "excluded"))

  # Add all thresholds into a single dataframe
  threshold_df <- data.frame(
    threshold = sapply(methods, function(m) RD_result$thresholds[[m]]$threshold),
    method = methods
  )
  # Create the plot (no faceting for single dataset)
  p <- ggplot(df_RD, aes(x = .data$distance)) +
    geom_histogram(aes(y = after_stat(density), fill = observation, color = observation),
                   binwidth = binwidth, alpha = 0.5) +
    geom_vline(
      data = threshold_df,
      aes(xintercept = threshold, color = method),
      linewidth = 0.8
    ) +

    scale_fill_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
    scale_color_manual(values = c(
      "SI" = "#1b9e77",
      "SI_boot" = "#2166ac",
      "MI" = "#7570b3",
      "MI_boot" = "#e7298a"
    )) +
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

#' Plot Robust Distance Histogram with Threshold
#'
#' Creates a histogram of robust distances with threshold line for outlier detection.
#'
#' @inheritParams thresh_result
#' @inheritParams RD_obj
#' @inheritParams alpha
#' @inheritParams binwidth
#' @inheritParams ...
#'
#' @return A ggplot object with histogram colored by inclusion status and threshold line.
#' @keywords internal
plot_RD_histogram <- function(thresh_result, RD_obj, alpha = 0.01, binwidth = 0.1, ...) {

  stopifnot("thresh_result must have a $threshold element" = !is.null(thresh_result$threshold))
  stopifnot("RD_obj must have $RD and $ind_incld elements" = !is.null(RD_obj$RD) && !is.null(RD_obj$ind_incld))

  # Extract what we need from the pre-computed objects
  RD <- RD_obj$RD
  ind_incld <- RD_obj$ind_incld
  t <- length(RD)

  # Get threshold and scaling from SI_result
  q_HR <- thresh_result$threshold
  scale_factor <- 1

  # Create observation status
  observation <- ifelse(seq_len(t) %in% ind_incld, "included", "excluded")

  # Extend tail range
  excl_vals <- RD[observation == "excluded"] * scale_factor
  dens_excld <- density(excl_vals, from = 0)
  peak_density <- max(dens_excld$y)
  xvals <- seq(0, max(c(excl_vals, q_HR), na.rm = TRUE) * 1.1, length.out = 512)

  # Create data frames for plotting (single dataset, no "data" column needed)
  df_RD <- data.frame(
    distance = RD * scale_factor,
    observation = observation
  )

  df_HR <- data.frame(x = xvals, y = rep(0, length(xvals)), q_HR = q_HR)
  df_HR_label <- data.frame(
    x = q_HR,
    y = peak_density * 1.05,
    label = paste0("")
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
    label = ""
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
    ggplot2::labs(title = "Outlier Matrix", x = "Variable", y = "Observations") +
    ggplot2::theme(aspect.ratio = 1.2) #+
    # scale_x_continuous(expand = c(0.01, 0)) +
    # scale_y_continuous(expand = c(0.01, 0))

  return(p)
}

#' Plot Multiple Imputation Results from RD Analysis
#'
#' Creates time series plots showing original data, temporal imputation,
#' and multiple imputation results with outlier locations highlighted.
#'
#' @param x An object of class "RD" from RD() or threshold_RD().
#' @param ... Additional arguments (unused).
#'
#' @return Prints ggplot objects for each variable showing imputation results.
#' @keywords internal
plot_imputations <- function(x) {
  outliers = TRUE # Want red dots
  x_data <- attr(x, "x_data")
  out_result <- attr(x, "univOut_hk")
  imp_result <- attr(x, "impute_univOut_hk")
  multiple_imp <- attr(x, "MImpute_hk")
  stopifnot("Multiple imputation data not available. This plot only available for MI and MI_boot methods." =
              !is.null(attr(x, "MImpute_hk")))
  NA_locs <- imp_result$NA_locs
  imputed_datasets <- multiple_imp$imp_datasets
  num_imputations <- length(imputed_datasets)
  num_columns <- ncol(out_result$outliers)  # Number of time series (columns)
  # # Set seed for reproducibility in selecting a random imputed dataset
  # set.seed(123)
  # Choose a random imputed dataset to highlight in orange
  random_imputed_idx <- sample(1:num_imputations, 1)
  method_used <- if ("interp" %in% as.character(x$call)) "interp" else "mean"
  message("Using method: ", method_used)
  # Loop through each column (IC)
  for (q in 1:num_columns) {
    # Prepare data for the current column (IC)
    plot_data <- data.frame(Time = 1:nrow(x_data), Original = x_data[, q])
    # Find the NA locations for this column (IC q)
    na_rows <- NA_locs[NA_locs[, 2] == q, 1]
    # Create the base plot without the original black line yet
    p <- ggplot()
    if (outliers) {
      # Add red dots for the missing data points (NA locations)
      p <- p + geom_point(data = data.frame(Time = na_rows, Original = x_data[na_rows, q]),
                          aes(x = Time, y = Original, color = "Outliers"), size = 2, alpha = 0.7)
    }
    # Add temporal imputation line (GREEN DASHED)
    p <- p + geom_line(data = data.frame(Time = 1:nrow(imp_result$imp_data), Temporal = imp_result$imp_data[, q]),
                       aes(x = Time, y = Temporal, color = "Temporal Imputation"), size = 1.0, linetype = "dashed")
    # Add imputed datasets as blue lines with transparency
    for (i in 1:num_imputations) {
      imp_data <- imputed_datasets[[i]]
      imp_data_plot <- data.frame(Time = 1:nrow(imp_data), Imputed = imp_data[, q])
      # If the current imputation is the randomly selected one, color it orange
      if (i == random_imputed_idx) {
        p <- p + geom_line(data = imp_data_plot, aes(x = Time, y = Imputed, color = "Random MI"),
                           size = 0.8)
      } else {
        p <- p + geom_line(data = imp_data_plot, aes(x = Time, y = Imputed, color = "Multiple Imputations"),
                           size = 0.5, alpha = 0.02)
      }
    }
    # Add the original time series in black (plotted last)
    if (method_used == "mean") {
      p <- p + geom_point(data = plot_data, aes(x = Time, y = Original, color = "Original"), size = 0.8)
    } else {
      p <- p + geom_line(data = plot_data, aes(x = Time, y = Original, color = "Original"), size = 0.8)
    }

    # Add manual color scale to match original colors
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
  }
}
