#' Plot Robust Distances with Custom Thresholds from Multiple Methods
#'
#' Plots robust distance histograms with color-coded inclusion status
#' and vertical lines for multiple thresholding methods (MI, MI_boot, etc.).
#'
#' @param results A named list where each element is a dataset (e.g., abide1, abide2),
#' containing precomputed RD and thresholds as in your `results$abide1$RD_org`, etc.
#' @return A ggplot object.
#' @export
plot_RD_with_thresholds <- function(results) {
  dataset_names <- names(results)
  
  # Step 1: Combine RD + labels
  df_combined <- do.call(rbind, lapply(dataset_names, function(nm) {
    df <- data.frame(
      RD = results[[nm]]$RD_org$RD,
      label = ifelse(seq_along(results[[nm]]$RD_org$RD) %in% results[[nm]]$RD_org$ind_incld, 
                     "included", "excluded"),
      dataset = toupper(nm)  # to match "ABIDE1", "ABIDE2"
    )
    return(df)
  }))
  
  df_combined$label <- factor(df_combined$label, levels = c("included", "excluded"))
  df_combined$dataset <- factor(df_combined$dataset, levels = toupper(dataset_names))
  
  # Step 2: Thresholds
  threshold_df <- data.frame(
    dataset = rep(toupper(dataset_names), each = 4),
    threshold = c(
      results[[1]]$MI$LB95_CI,
      results[[1]]$MI_boot$final_threshold,
      results[[1]]$SI$SI_threshold,
      results[[1]]$SI_boot$LB_CI,
      results[[2]]$MI$LB95_CI,
      results[[2]]$MI_boot$final_threshold,
      results[[2]]$SI$SI_threshold,
      results[[2]]$SI_boot$LB_CI
    ),
    method = rep(c("MI", "MI_boot", "SI", "SI_boot"), times = length(results))
  )
  
  # Step 3: Final plot
  p <- ggplot(df_combined, aes(x = RD, fill = label)) +
    geom_histogram(aes(y = ..density..), bins = 30, position = "identity", color = "black") +
    geom_vline(data = threshold_df, aes(xintercept = threshold, linetype = method, color = method),
               linewidth = 0.8, show.legend = TRUE) +
    facet_wrap(~ dataset, scales = "free") +
    scale_x_log10() +
    scale_fill_manual(values = c("included" = "#66c2a5", "excluded" = "#fc8d62")) +
    scale_linetype_manual(
      name = "Threshold",
      values = c("MI" = "solid", "MI_boot" = "dashed", "SI" = "dotdash", "SI_boot" = "twodash")
    ) +
    scale_color_manual(
      name = "Threshold",
      values = c("MI" = "blue", "MI_boot" = "magenta", "SI" = "green4", "SI_boot" = "darkorange")
    ) +
    labs(
      x = "Robust Distance (log10)",
      y = "Density",
      fill = "Observations"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      strip.text = element_text(face = "bold", size = 14)
    )
  
  return(p)
}
