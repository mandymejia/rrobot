#' Plot Robust Mahalanobis distances and theoretical threshold from F-distribution
#'
#' @param data A numeric matrix (T × Q) containing multivariate observations (rows = time/obs, cols = variables).
#' @param alpha Numeric in (0, 1). Quantile level used to compute the threshold via tresh_F(). Default 0.01.
#' @param binwidth Bin width for the histogram of robust distances. Default 0.1.
#' @param log Logical. If TRUE, applies log10 scale to the x-axis (on 1 + x). Default FALSE.
#' @param show_f_density Logical. Draw scaled F-density curve. Default TRUE.
#' @return A ggplot object.
#' @seealso comp_RD, tresh_F
#' @export
plot_RD <- function(data,
                    alpha = 0.01,
                    binwidth = 0.1,
                    log = FALSE,
                    show_f_density = TRUE) {

  # --- checks ---
  if (!is.matrix(data) || !is.numeric(data))
    stop("`data` must be a numeric matrix (T × Q).")
  if (any(!is.finite(colSums(is.na(data)))))
    warning("Data contain NA values; comp_RD() must handle them appropriately.")

  Tn <- nrow(data)
  Q  <- ncol(data)

  # --- robust distances & inclusion mask ---
  rob_obj    <- comp_RD(data)
  RD         <- rob_obj$RD
  ind_incld  <- rob_obj$ind_incld
  observation <- ifelse(seq_len(Tn) %in% ind_incld, "included", "excluded")

  # --- F-threshold via tresh_F ---
  fitF_obj <- tresh_F(Q, n = Tn, h = rob_obj$h, quantile = alpha)
  q_HR     <- fitF_obj$threshold * fitF_obj$scale
  scale_x  <- fitF_obj$scale
  df1      <- fitF_obj$df[1]
  df2      <- fitF_obj$df[2]

  # --- density support (use excluded tail + threshold) ---
  excl_vals <- (RD[observation == "excluded"]) * scale_x
  if (length(excl_vals) < 2L || !any(is.finite(excl_vals))) {
    # Fallback: use all RD if excluded is empty/degenerate
    excl_vals <- RD * scale_x
  }
  dens_excld <- density(excl_vals[is.finite(excl_vals)], from = 0)
  peak_y     <- max(dens_excld$y)
  x_max      <- max(c(excl_vals, q_HR), na.rm = TRUE) * 1.1
  xvals      <- seq(0, x_max, length.out = 512)

  # --- scaled F curve to match empirical peak ---
  f_y_raw <- df(xvals, df1 = df1, df2 = df2)
  scale_y <- peak_y / max(f_y_raw)
  yvals   <- f_y_raw * scale_y

  # --- assemble frames for plotting ---
  df_hist <- data.frame(
    distance    = RD * scale_x,
    observation = factor(observation, levels = c("included", "excluded"))
  )
  df_fcurve <- data.frame(x = xvals, y = yvals, q_HR = q_HR)
  df_label  <- data.frame(
    x = q_HR,
    y = max(yvals, na.rm = TRUE) * 1.05,
    label = paste0("F(", 1 - alpha, ")")
  )

  # --- ggplot ---
  library(ggplot2)
  p <- ggplot(df_hist, aes(x = distance)) +
    geom_histogram(aes(y = after_stat(density), fill = observation, color = observation),
                   binwidth = binwidth, alpha = 0.5) +
    geom_vline(xintercept = q_HR, linetype = "dashed", linewidth = 0.8, color = "#333333") +
    { if (show_f_density) geom_line(data = df_fcurve, aes(x = x, y = y),
                                    inherit.aes = FALSE, linewidth = 0.8, color = "#333333") } +
    geom_text(data = df_label, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = -0.1, vjust = -0.5, size = 4.5, color = "#333333") +
    scale_fill_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
    scale_color_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
    { if (log) scale_x_log10() else scale_x_continuous(n.breaks = 10) } +
    labs(y = "Density",
         x = if (log) "log10(Squared robust distances)" else "Squared robust distances") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      strip.text = element_text(size = 13),
      legend.key = element_blank()
    )

  return(p)
}
