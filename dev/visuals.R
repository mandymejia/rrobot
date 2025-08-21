# visualizations of threshold values using the rrobot outputs directly
# lets load the object from the local folder
library(here)
abide1_obj <- readRDS(here::here("dev", "fixtures", "ABIDE1.rds"))
RD_abide1 <- abide1_obj$RD_obj$RD
RD_sh_abide1 <- abide1_obj$RD_obj_shash$RD
xmax <- log(max(RD_sh_abide1))
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
library(ggplot2)
source(here::here("dev", "ICA_extract_kurt.R"))
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

df_qq <- rbind(
  data.frame(Theoretical = qq_before$x,
             Sample = qq_before$y,
             Type = "Before"),
  data.frame(Theoretical = qq_after$x,
             Sample = qq_after$y,
             Type = "After")
)

# Symmetric axis limits
lim <- max(abs(c(df_qq$Theoretical, df_qq$Sample)))

ggplot(df_qq, aes(x = Theoretical, y = Sample, color = Type)) +
  geom_point(alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  coord_equal(xlim = c(-3,3), ylim = c(-3, 3)) +
  labs(title = "QQ Plot: Before vs After SHASH Transformation",
       x = "Theoretical Quantiles (N(0,1))",
       y = "Sample Quantiles") +
  scale_color_manual(values = c("Before" = "#FF99FF",
                                "After" = "lightgreen")) +
  theme_minimal(base_size = 14)



library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# Ensure columns have names
if (is.null(colnames(hk))) colnames(hk) <- paste0("V", seq_len(ncol(hk)))
if (is.null(colnames(Zhk))) colnames(Zhk) <- colnames(hk)

# Helper: make QQ df for one vector
qq_df_one <- function(x, type_label) {
  q <- qqnorm(x, plot.it = FALSE)
  data.frame(Theoretical = q$x, Sample = q$y, Type = type_label, stringsAsFactors = FALSE)
}

scale_MAD <- function(x){
  x_med <- median(x, na.rm = TRUE)
  mad_val <- 1.4826 * median(abs(x - x_med), na.rm = TRUE)
  z <- (x - x_med) / mad_val
}


# Build combined QQ data for all columns
df_qq_all <- map_dfr(colnames(hk), function(nm){
  before <- qq_df_one(scale_MAD(hk[, nm]),  "Before(ABIDE1)")
  after  <- qq_df_one(Zhk[, nm], "After (SHASH-Normal)")
  bind_rows(before, after) %>% mutate(Column = nm)
})

# to order variables from 1 to 25
df_qq_all$Column <- factor(
  df_qq_all$Column,
  levels = paste0("V", seq_len(ncol(hk)))  # forces V1, V2, ..., V25 order
)

# Symmetric limits (global); adjust 'cap' if you prefer wider/tighter view
cap <- 4
lim <- cap

ggplot(df_qq_all, aes(x = Theoretical, y = Sample, color = Type)) +
  geom_point(alpha = 0.6, size = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.6) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  facet_wrap(~ Column, scales = "fixed") +
  labs(title = "QQ Plots: Before vs After SHASH Transformation (All Columns)",
       x = "Theoretical Quantiles (N(0,1))",
       y = "Sample Quantiles") +
  scale_color_manual(values = c(
    "Before(ABIDE1)"   = "#990000",
    "After (SHASH-Normal)"  = "#339900"
  )) +
  theme_minimal(base_size = 13)


###  adding the F density line over the RD histogram for SHASH
# truncate the values
Zhk[ Zhk > 100] <- 100
Zhk[ Zhk < -100] <- -100

rrobot::plot_RD(data_list = Zhk, log = TRUE)
FP_plotRD <- function(data,
                      alpha = 0.01,
                      binwidth = 0.1,
                      log = FALSE,
                      show_f_density = TRUE) {

  # --- checks ---
  if (!is.matrix(data) || !is.numeric(data))
    stop("`data` must be a numeric matrix (T × Q).")
  if (any(!is.finite(colSums(is.na(data)))))
    warning("Data contain NA values; compute_RD() must handle them appropriately.")

  Tn <- nrow(data)
  Q  <- ncol(data)

  # --- robust distances & inclusion mask ---
  rob_obj    <- compute_RD(data)
  RD         <- rob_obj$RD
  ind_incld  <- rob_obj$ind_incld
  observation <- ifelse(seq_len(Tn) %in% ind_incld, "included", "excluded")

  # --- F-threshold via thresh_F ---
  fitF_obj <- thresh_F(Q, n = Tn, h = rob_obj$h, quantile = alpha)
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

FP_plotRD(Zhk, log = TRUE)

