# visualizations of threshold values using the rrobot outputs directly
# lets load the object from the local folder
abide1_obj <- readRDS("~/Documents/GitHub/RobOutlier-paper/results/ABIDE/ABIDE1.rds")

RD_abide1 <- abide1_obj$RD_obj$RD
RD_sh_abide1 <- abide1_obj$RD_obj_shash$RD

df1 <- data.frame(
  RD = abide1_obj$RD_obj$RD,
  label = ifelse(seq_along(abide1_obj$RD_obj$RD) %in% abide1_obj$RD_obj$ind_incld, "included", "excluded")
)

df1$label <- factor(df1$label, levels = c("included", "excluded"))

df_thresh <- data.frame(
  threshold = c(
    abide1_obj$thresholds$SI$SI_threshold,
    abide1_obj$thresholds$SI_boot$LB_CI,
    abide1_obj$thresholds$MI$LB95_CI,
    abide1_obj$thresholds$MI_boot$final_threshold,
    abide1_obj$thresholds$F$threshold),
  method = c("SI", "SI_boot", "MI", "MI_boot", "F")
)


# Final ggplot
library(ggplot2)
ggplot(df1, aes(x = RD, fill = label)) +
  geom_histogram(aes(y = ..density..), bins = 40, position = "identity", alpha= 0.9) +
  geom_vline(
    data = df_thresh,
    aes(xintercept = threshold, color = method),
    linewidth = 0.8
  ) +
  scale_x_log10() +
  scale_fill_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
  scale_color_manual(
    name = "Threshold",
    values = c("SI" = "#1b9e77", "SI_boot" = "#d95f02", "MI" = "#7570b3", "MI_boot" = "#e7298a", "F" = "red")
  ) +
  geom_text(
    data = df_thresh,
    aes(x = threshold, y = Inf, label = method, color = method),
    angle = 90, vjust = -0.9, hjust = 5, size = 4,
    inherit.aes = FALSE, show.legend = FALSE
  )+
  labs(
    x = "Robust Distance",
    y = "Density",
    fill = "Observations"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text = element_text(face = "bold", size = 14)
  )

