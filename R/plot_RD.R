#' Plot Robust Distances and HR Thresholds for One or Multiple Datasets
#'
#' Generates histogram(s) of squared robust Mahalanobis distances computed from robust covariance estimates.
#' Overlays vertical threshold lines (e.g., HR threshold) and, optionally, F-distribution curves derived from Hotelling’s T² theory.
#' Works with either a single dataset or a list of datasets.
#'
#' @param data_list A numeric matrix (T × Q), or a named list of such matrices. Each matrix should contain multivariate observations.
#' @param alpha Numeric in (0, 1). The quantile level used to compute the HR threshold via \code{\link{outlier_F}} if \code{HR} is not provided. Default is 0.01.
#' @param HR Optional numeric vector of threshold(s) for robust distances. If \code{NULL} (default), the threshold is computed via \code{\link{outlier_F}} for each dataset. If a single value is provided, it will be applied to all datasets. If a named vector is provided, names must match dataset names.
#' @param n_rep Integer. Number of repetitions used to estimate the mean false positive rate (FPR) by recomputing the RD and comparing to the threshold. Default is 1000.
#' @param label_map Optional named character vector to override facet labels. For example, \code{c("Gaussian" = "N(0,1)", "Skewed" = "Gamma(5,1)")} will relabel facet strips.
#' @param binwidth Bin width for the histogram of robust distances. Default is 0.1.
#'
#' @return A \code{ggplot} object that includes:
#' \itemize{
#'   \item Histogram of robust distances by inclusion status (included vs. excluded)
#'   \item Optional F-distribution density curve (if \code{HR} is \code{NULL})
#'   \item Vertical dashed line for the supplied or computed threshold
#'   \item FPR annotation per dataset panel
#' }
#'
#' @seealso \code{\link{comp_RD}}, \code{\link{outlier_F}}
#' @importFrom ggplot2 ggplot aes geom_histogram geom_line geom_vline geom_text
#'   scale_fill_manual scale_color_manual coord_cartesian theme_minimal
#'   facet_grid after_stat theme element_text element_blank
#'   scale_x_continuous as_labeller label_value
#' @importFrom stats df
#' @export
plot_RD <- function(data_list
                    , alpha = 0.01
                    , HR = NULL
                    , n_rep = 1000
                    , label_map = NULL
                    , binwidth = 0.1) {
  if (is.matrix(data_list)) {
    data_list <- list(Data = data_list)
  }
  if (is.null(names(data_list))) {
    names(data_list) <- paste0("Data", seq_along(data_list))
  }

  datanames <- names(data_list)
  t <- nrow(data_list[[1]])
  Q <- ncol(data_list[[1]])

  df_RD <- df_HR <- mean_FPs <- NULL

  for (i in seq_along(datanames)) {
    data <- data_list[[i]]
    dist_name <- datanames[i]
    message("Processing: ", dist_name)

    rob_obj <- comp_RD(data)
    RD <- rob_obj$RD
    ind_incld <- rob_obj$ind_incld
    ind_excld <- rob_obj$ind_excld
    observation <- ifelse(seq_len(t) %in% ind_incld, "included", "excluded")

    # Get threshold
    if (is.null(HR)) {
      fitF_obj <- outlier_F(Q, n = t, h = rob_obj$h, quantile = alpha)
      q_HR <- fitF_obj$threshold * fitF_obj$scale
      scale_factor <- fitF_obj$scale
      xvals <- seq(0, max(RD[ind_excld] * scale_factor, na.rm = TRUE), length.out = t)
      yvals <- df(xvals, df1 = fitF_obj$df[1], df2 = fitF_obj$df[2])
    } else {
      if (length(HR) == 1) {
        q_HR <- HR
      } else {
        if (is.null(names(HR)) || !(dist_name %in% names(HR))) {
          stop("If HR is a named vector, it must contain names matching each dataset.")
        }
        q_HR <- HR[dist_name]
      }
      scale_factor <- 1
      xvals <- seq(0, max(RD[ind_excld], na.rm = TRUE), length.out = t)
      yvals <- rep(NA, length(xvals))
    }

    df_HR <- rbind(df_HR, data.frame(
      x = xvals,
      y = yvals * 2,
      q_HR = q_HR,
      data = dist_name
    ))

    df_RD <- rbind(df_RD, data.frame(
      distance = RD * scale_factor,
      observation = observation,
      data = dist_name
    ))

    # Estimate false positives
    FPs_i <- replicate(n_rep, {
      rob_tmp <- comp_RD(data)
      thresh <- if (is.null(HR)) {
        outlier_F(Q, t, rob_tmp$h, quantile = alpha)$threshold
      } else if (length(HR) == 1) {
        HR
      } else {
        HR[dist_name]
      }
      mean(rob_tmp$RD >= thresh)
    })
    mean_FPs <- rbind(mean_FPs, data.frame(FPs = round(mean(FPs_i), 3), data = dist_name))
  }

  df_RD_HR <- merge(df_RD, df_HR[, c("x", "y", "q_HR", "data")], by = "data")
  df_RD_HR$observation <- factor(df_RD_HR$observation, levels = c("included", "excluded"))
  df_RD_HR$data <- factor(df_RD_HR$data, levels = datanames)

  anno_df <- data.frame(data = datanames, label = mean_FPs$FPs)
  facet_labeller <- if (!is.null(label_map)) as_labeller(label_map) else label_value

  p <- ggplot(df_RD_HR, aes(x = distance)) +
    geom_histogram(aes(y = after_stat(density), fill = observation, color = observation),
                   binwidth = binwidth, alpha = 0.5) +
    geom_vline(aes(xintercept = q_HR), linetype = "dashed", linewidth = 0.8, color = '#333333') +
    geom_text(data = anno_df, aes(label = label, x = df_HR$q_HR, y = Inf),
              hjust = 1.1, vjust = 1.1, color = "#333333",
              inherit.aes = FALSE, size = 5)

  if (is.null(HR)) {
    p <- p + geom_line(aes(x = x, y = y), size = 0.8, color = '#333333')
  }

  p <- p +
    facet_grid(data ~ ., labeller = facet_labeller, scales = "free_y") +
    scale_fill_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
    scale_color_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
    coord_cartesian(ylim = c(0, 2), xlim = c(0, 10)) +
    scale_x_continuous(n.breaks = 10) +
    ylab("Density") + xlab("Squared robust distances") +
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
