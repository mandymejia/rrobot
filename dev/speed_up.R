library(fMRIscrub)
data_matrix <- fMRIscrub::Dat1
# source("point to scripts read the funct)"
kurt_data <- ICA_extract_kurt(time_series = data_matrix) # use fMRIscrub
hk_abide1 <- kurt_data$hk
lk_abide1 <- kurt_data$lk
dim(hk_abide1)

# want to test if we need to save RDs for each methods:

t <- system.time({
  hk1_obj <- rrobot::compute_RD(hk_abide1, mode = "auto")
  hk1_thrs <- rrobot::threshold_RD(
                  x = hk_abide1, w = lk_abide1,
                  method = "all", RD_obj = hk1_obj,
                  impute_method = "interp",
                  M = 100, k = 5, B = 1000)
})

print(t)              # user / system / elapsed
t["elapsed"]          # numeric seconds

str(hk1_thrs)
hist(log(hk1_thrs$RD_obj$RD), breaks = 30)


library(ggplot2)

# --- Data for histogram -------------------------------------------------------
RD      <- hk1_thrs$RD_obj$RD
ind_in  <- hk1_thrs$RD_obj$ind_incld
label   <- rep("excluded", length(RD))
label[ind_in] <- "included"

eps <- 1e-12
z   <- log(pmax(RD, eps))
df  <- data.frame(z = z, label = factor(label, levels = c("included","excluded")))

# --- Collect thresholds (on RD scale), EXCLUDING SHASH ------------------------
thr_SI      <- hk1_thrs$thresholds$SI$threshold
thr_SI_boot <- hk1_thrs$thresholds$SI_boot$threshold
thr_MI      <- hk1_thrs$thresholds$MI$threshold
thr_MI_boot <- hk1_thrs$thresholds$MI_boot$threshold
# F threshold needs scaling to RD units:
thr_F       <- hk1_thrs$thresholds$F$threshold * hk1_thrs$thresholds$F$scale

df_thresh <- data.frame(
  method  = factor(c("SI (99%)", "SI_boot (2.5% CI LB)", "MI (2.5%)", "MI_boot", "F"),
                   levels = c("SI (99%)", "SI_boot (2.5% CI LB)", "MI (2.5%)", "MI_boot", "F")),
  thr_log = log(pmax(c(thr_SI, thr_SI_boot, thr_MI, thr_MI_boot, thr_F), eps))
)

# --- Plot ---------------------------------------------------------------------
ggplot(df, aes(x = z, fill = label)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, position = "identity", alpha = 0.6, color = "black") +
  geom_vline(data = df_thresh,
             aes(xintercept = thr_log, color = method),
             linewidth = 0.9, linetype = "solid") +
  scale_fill_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
  labs(
    title = "Histogram of log(RD) with Thresholds (non-SHASH)",
    x = "log(RD)",
    y = "Density",
    fill = "Observation",
    color = "Threshold"
  ) +
  theme_minimal(base_size = 12)
