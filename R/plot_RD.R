#' Plot Robust Mahalanobis distances and theoretical threshold from F-distribution
#'
#' Generates histogram(s) of squared robust Mahalanobis distances from robust covariance estimates.
#' Overlays vertical threshold lines (computed via \code{\link{thresh_F}}) and a scaled F-distribution curve.
#' Optionally applies log-scaling to the x-axis for better visualization of heavy-tailed data.
#'
#' @param data_list A numeric matrix (T × Q), or a named list of such matrices. Each matrix should contain multivariate observations.
#' @param alpha Numeric in (0, 1). The quantile level used to compute the threshold via \code{\link{thresh_F}}. Default is 0.01.
#' @param binwidth Bin width for the histogram of robust distances. Default is 0.1.
#' @param log Logical. If TRUE, applies log1p() transformation to the x-axis. Default is FALSE.
#' @param show_f_density Logical. If \code{TRUE}, overlay a scaled F-distribution density curve. Default \code{TRUE}.
#'
#' @return A \code{ggplot} object that includes:
#' \itemize{
#'   \item Histogram of robust distances by inclusion status (included vs. excluded)
#'   \item Scaled F-distribution density curve
#'   \item Vertical dashed line for the computed threshold
#'   \item Optional F(α) annotation next to the threshold
#' }
#'
#' @seealso \code{\link{RD}}, \code{\link{thresh_F}}
#' @importFrom ggplot2 ggplot aes geom_histogram geom_line geom_vline geom_text
#'   scale_fill_manual scale_color_manual theme_minimal facet_grid after_stat theme element_text
#'   element_blank scale_x_continuous scale_x_log10 ylab xlab
#' @importFrom stats df density
#' @importFrom rlang .data
#' @export
plot_RD <- function(data_list,
                    alpha = 0.01,
                    binwidth = 0.1,
                    log = FALSE,
                    show_f_density = TRUE) {
  if (is.matrix(data_list)) {
    data_list <- list(Data = data_list)
  }
  if (is.null(names(data_list))) {
    names(data_list) <- paste0("Data", seq_along(data_list))
  }

  datanames <- names(data_list)
  t <- nrow(data_list[[1]])
  Q <- ncol(data_list[[1]])

  df_RD <- df_HR <- df_HR_label <- NULL

  for (i in seq_along(datanames)) {
    data <- data_list[[i]]
    dist_name <- datanames[i]
    message("Processing: ", dist_name)

    rob_obj <- RD(data)
    RD <- rob_obj$RD
    ind_incld <- rob_obj$ind_incld
    observation <- ifelse(seq_len(t) %in% ind_incld, "included", "excluded")

    # Compute threshold from thresh_F
    fitF_obj <- thresh_F(Q, n = t, h = rob_obj$h, quantile = alpha)
    q_HR <- fitF_obj$threshold * fitF_obj$scale
    scale_factor <- fitF_obj$scale
    df1 <- fitF_obj$df[1]
    df2 <- fitF_obj$df[2]

    # Extend tail range
    excl_vals <- RD[observation == "excluded"] * scale_factor
    dens_excld <- density(excl_vals, from = 0)
    peak_density <- max(dens_excld$y)
    # x_max <- max(quantile(excl_vals, 0.99, na.rm = TRUE), q_HR) * 2.5
    # xvals <- seq(0, x_max, length.out = 512)
    xvals <- seq(0, max(c(excl_vals, q_HR), na.rm = TRUE) * 1.1, length.out = 512)

    # F-distribution curve scaled to match empirical peak
    f_density <- df(xvals, df1 = df1, df2 = df2)
    scale_factor_y <- peak_density / max(f_density)
    yvals <- f_density * scale_factor_y

    # Store all data
    df_HR <- rbind(df_HR, data.frame(x = xvals, y = yvals, q_HR = q_HR, data = dist_name))
    df_HR_label <- rbind(df_HR_label, data.frame(
      x = q_HR,
      y = max(yvals, na.rm = TRUE) * 1.05,
      label = paste0("F(", 1 - alpha, ")"),
      data = dist_name
    ))

    df_RD <- rbind(df_RD, data.frame(
      distance = RD * scale_factor,
      observation = observation,
      data = dist_name
    ))
  }

  # Merge for plotting
  df_RD_HR <- merge(df_RD, df_HR[, c("x", "y", "q_HR", "data")], by = "data")
  df_RD_HR$observation <- factor(df_RD_HR$observation, levels = c("included", "excluded"))
  df_RD_HR$data <- factor(df_RD_HR$data, levels = datanames)
  df_HR_label$data <- factor(df_HR_label$data, levels = datanames)

  # Conditional faceting
  facet_layer <- if (length(data_list) > 1) {
    facet_grid(data ~ ., scales = "free_y")
  } else NULL

  df_RD_HR <- df_RD_HR[df_RD_HR$distance > 0, ]
  df_HR <- df_HR[df_HR$q_HR > 0, ]
  df_HR_label <- df_HR_label[df_HR_label$x > 0, ]

  # Plot
  p <- ggplot(df_RD_HR, aes(x = .data$distance)) +
    geom_histogram(aes(y = after_stat(density), fill = observation, color = observation),
                   binwidth = binwidth, alpha = 0.5) +
    geom_vline(aes(xintercept = q_HR), linetype = "dashed", linewidth = 0.8, color = '#333333') +
    (if (show_f_density) geom_line(aes(x = .data$x, y = .data$y), size = 0.8, color = '#333333') else NULL) +
    geom_text(data = df_HR_label,
              aes(x = .data$x, y = .data$y, label = .data$label),
              inherit.aes = FALSE,
              hjust = -0.1, vjust = -0.5, size = 4.5, color = "#333333") +
    facet_layer +
    # coord_cartesian(xlim = c(NA, max_display)) +
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
