abide1shash_obj <- readRDS("~/Documents/GitHub/RobOutlier-paper/results/ABIDE/ABIDE1_SHASH.rds")
df_RD <- readRDS("~/Documents/GitHub/RobOutlier-paper/results/df_RD.rds")
df_HR <- readRDS("~/Documents/GitHub/RobOutlier-paper/results/df_HR.rds")

data <- abide1shash_obj$x_norm
Tn <- nrow(data)
Q  <- ncol(data)

# --- robust distances & inclusion mask ---
RD         <- abide1shash_obj$RD_obj_shash$RD
ind_incld  <- abide1shash_obj$RD_obj_shash$ind_incld
observation <- ifelse(seq_len(Tn) %in% ind_incld, "included", "excluded")

# --- F-threshold via Fit_F ---
q_HR     <- abide1shash_obj$thresholds$threshold
scale_x  <- abide1shash_obj$thresholds$scale
df1      <- abide1shash_obj$thresholds$df[1]
df2      <- abide1shash_obj$thresholds$df[2]


RD_out <- RD[ind_excld_shash]
xvals <- seq(0, max(RD_out * scale_x), length.out = Tn)
yvals <- df(xvals, df1 = df[1], df2 = df[2])

# --- assemble frames for plotting ---
df_hist <- data.frame(
  distance    = RD * scale_x,
  observation = factor(observation, levels = c("included", "excluded"))
)
df_fcurve <- data.frame(x = xvals, y = yvals, q_HR = q_HR*scale_x)

# Find max histogram density
max_hist_density <- max(
  ggplot_build(
    ggplot(df_hist, aes(x = distance)) +
      geom_histogram(aes(y = after_stat(density)), binwidth = 0.3)
  )$data[[1]]$density
)

# Scale your fcurve's y to match histogram density scale
df_fcurve <- df_fcurve %>%
  mutate(y_scaled = y * (max_hist_density / max(y, na.rm = TRUE)))

# --- ggplot ---
library(ggplot2)
ggplot(df_hist, aes(x = distance))+
  geom_histogram(aes(y = after_stat(density), fill = observation, color = observation),
                 binwidth = 0.3, alpha = 0.5) +
  geom_vline(xintercept = q_HR, linetype = "dashed", linewidth = 0.8, color = "#333333") +
  scale_x_log10() +
  geom_line(data = df_fcurve, aes(x = x, y = y_scaled),
                                  inherit.aes = FALSE, linewidth = 0.8, color = "#333333") +
  scale_fill_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
  scale_color_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
  labs(y = "Density",
       x = "log10(Squared robust distances") +
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




## --- Append ABIDE1 (SHASH) to the existing three panels -----------------
abide1shash_obj <- readRDS("~/Documents/GitHub/RobOutlier-paper/results/ABIDE/ABIDE1_SHASH.rds")

# data + RD on original scale used in your first panels
RD_abide   <- abide1shash_obj$RD_obj_shash$RD
ind_incld  <- abide1shash_obj$RD_obj_shash$ind_incld
Tn         <- length(RD_abide)

obs_abide <- rep("excluded", Tn)
obs_abide[ind_incld] <- "included"

# F params from object
df1 <- abide1shash_obj$thresholds$df[1]
df2 <- abide1shash_obj$thresholds$df[2]
s   <- abide1shash_obj$thresholds$scale
thr <- abide1shash_obj$thresholds$threshold * s

RD_out_abide <- RD_abide[obs_abide == "excluded"]
xmax_abide   <- max(RD_out_abide, na.rm = TRUE)
log_scale <- log(seq(0.01, xmax_abide, length.out = 1000))
xvals_abide  <- log_scale
yvals_abide  <- (1/s) * stats::df(xvals_abide, df1 = df1, df2 = df2)

# To match your earlier visual scaling (you multiplied by 2):
yvals_abide_scaled <-  yvals_abide

# Frames matching your existing structure
tmp_HR_abide <- data.frame(
  x    = xvals_abide,
  y    = yvals_abide_scaled,
  q_HR = thr,
  data = "ABIDE1 (SHASH)"
)

tmp_RD_abide <- data.frame(
  distance    = RD_abide * s,
  observation = obs_abide,
  data        = "ABIDE1 (SHASH)"
)

# Bind to existing frames then rebuild df_RD_HR
df_HR <- rbind(tmp_HR_abide)
df_RD <- rbind(tmp_RD_abide)

df_RD_HR <- merge(df_RD, df_HR[, c("x", "y", "q_HR", "data")], by = "data")
df_RD_HR$observation <- factor(df_RD_HR$observation, levels = c("included", "excluded"))

library(ggplot2)

mean_emp <- mean(df_RD_HR$distance, na.rm = TRUE)

# Build histogram once to get its max density
p_tmp <- ggplot(df_RD_HR, aes(x = distance)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3)
max_hist_density <- max(ggplot_build(p_tmp)$data[[1]]$density, na.rm = TRUE)

# Scale the F curve up to match
curve_scale <- max_hist_density / max(tmp_HR_abide$y, na.rm = TRUE)
df_HR$y_scaled <- df_HR$y * curve_scale

# Plot using y_scaled
ggplot(df_RD_HR, aes(x = log(distance))) +
  geom_histogram(aes(y = after_stat(density), fill = observation, color = observation),
                 binwidth = 0.3, alpha = 0.5) +
  geom_line(data = df_HR, aes(x = x, y =y),
            inherit.aes = FALSE, linewidth = 0.8, color = "#333333") +
  geom_vline(aes(xintercept = q_HR), linetype = "dashed", linewidth = 0.8, color = "#333333")



