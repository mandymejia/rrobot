# visualizations of threshold values using the rrobot outputs directly
# lets load the object from the local folder
abide1_obj <- readRDS("~/Documents/GitHub/RobOutlier-paper/results/ABIDE/ABIDE1.rds")

RD_abide1 <- abide1_obj$RD_obj$RD
RD_sh_abide1 <- abide1_obj$RD_obj_shash$RD

df1 <- data.frame(
  RD = abide1_obj$RD_obj$RD,
  label = ifelse(seq_along(abide1_obj$RD_obj$RD) %in% abide1_obj$RD_obj$ind_incld, "included", "excluded")
)

df_comb <- rbind(df1, df2)
df_comb$label <- factor(df_comb$label, levels = c("included", "excluded"))


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
p1 <- ggplot(df1, aes(x = log(RD), fill = label)) +
  geom_histogram(aes(y = ..density..), bins = 40, position = "identity", alpha= 0.7, color = "black") +
  geom_vline(
    data = df_thresh,
    aes(xintercept = log(threshold), linetype = method, color = method),
    linewidth = 0.8
  ) +
  # scale_x_log10() +
  xlim(0, xmax) +
  scale_fill_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
  scale_color_manual(
    name = "Threshold",
    values = c("SI" = "black", "SI_boot" = "black", "MI" = "gray", "MI_boot" = "gray", "F" = "red")
  ) +
  scale_linetype_manual(
    name = "Threshold",
    values = c("SI" = 1, "SI_boot" = 2, "MI" = 1, "MI_boot" = 2, "F" = 1)
  ) +
  labs(
    x = "log of Robust Distance",
    y = "Density",
    fill = "Observations"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text = element_text(face = "bold", size = 14)
  )


##----- shash-------------------------------------------------------------------
RD_sh_abide1 <- abide1_obj$RD_obj_shash$RD

df2 <- data.frame(
  RD = RD_sh_abide1,
  label = ifelse(seq_along(RD_sh_abide1) %in% abide1_obj$RD_obj_shash$ind_incld,
                 "included", "excluded")
)

xmax <- log(max(RD_sh_abide1))

df2$label <- factor(df2$label, levels = c("included", "excluded"))
options(scipen = 999)
p2 <- ggplot(df2, aes(x = log(RD), fill = label)) +
  geom_histogram(aes(y = ..density..), bins = 40, position = "identity",
                 alpha = 0.7, color = "black") +
  geom_vline(xintercept = log(abide1_obj$thresholds$SHASH$threshold),
             color = "red", linewidth = 0.8, show.legend = FALSE) +
  # scale_x_log10() +
  xlim(0, xmax) +
  scale_fill_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
  labs(x = "log Robust Distance", y = "Density", fill = "Observations") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        strip.text = element_text(face = "bold", size = 14))

####------scatterplot-----------------------------------------------------------
library(dplyr)
library(ggplot2)

# Vectors of distances (same length n)
RD       <- abide1_obj$RD_obj$RD
RD_shash <- abide1_obj$RD_obj_shash$RD

id <- seq_along(RD)

stopifnot(length(RD) == length(RD_shash))

incl_RD  <- id %in% abide1_obj$RD_obj$ind_incld
incl_SH  <- id %in% abide1_obj$RD_obj_shash$ind_incld

df <- tibble(
  id, RD, RD_shash,
  incl_RD  = incl_RD,
  incl_SH  = incl_SH,
  status = case_when(
    incl_RD  &  incl_SH  ~ "Included by both",
    !incl_RD  & !incl_SH  ~ "Excluded by both",
    incl_RD  & !incl_SH  ~ "Only RD included",
    !incl_RD  &  incl_SH  ~ "Only SHASH included"
  )
) %>%
  # for log scales, remove non-positive or non-finite values
  filter(is.finite(RD), is.finite(RD_shash), RD > 0, RD_shash > 0) %>%
  mutate(status = factor(
    status,
    levels = c("Included by both", "Excluded by both", "Only RD included", "Only SHASH included")
  ))

ggplot(df, aes(x = RD, y = RD_shash, color = status)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +     # y = x reference
  geom_hline(yintercept = 100) +
  geom_vline(xintercept = 100) +
  geom_point(alpha = 0.8, size = 2) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c(
    "Included by both"  = "#1b9e77",
    "Excluded by both"  = "#808080",
    "Only RD included"  = "#d95f02",
    "Only SHASH included" = "#7570b3"
  )) +
  labs(
    x = "Robust Distance (RD)",
    y = "Robust Distance (SHASH)",
    color = "Agreement status",
    title = "RD vs RD_SHASH"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")


table(df$status)


gridExtra::grid.arrange(p1, p2, nrow= 2)

#################################################################################
## QQ plots before and after the SHASH transformation
#################################################################################
library(fMRIscrub)
data_matrix <- fMRIscrub::Dat1
kurt_data <- ICA_extract_kurt(time_series = data_matrix)
hk <- kurt_data$hk
dim(kurt_data$hk)
Zhk <- hk
Zhk[,] <- NA
for(i in seq_len(ncol(hk))){
  hk_i <- hk[,i]
  Zhk[,i] <- SHASH_out(hk_i)$x_norm
}

qq_before <- qqnorm(hk[,1], plot.it = FALSE)
qq_after  <- qqnorm(Zhk[,1], plot.it = FALSE)


