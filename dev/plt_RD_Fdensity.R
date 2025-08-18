library(ggplot2)
library(here)
# --- Load ---------------------------------------------------------------------
abide1shash_obj <- readRDS(here::here("dev", "fixtures", "ABIDE1.rds"))

RD_obj <- abide1shash_obj$RD_obj_shash
rd2    <- as.numeric(RD_obj$RD)          # RDsquare vector

Fpar <- abide1shash_obj$thresholds
df1  <- Fpar$df[1]
df2  <- Fpar$df[2]
s    <- 1 / Fpar$scale
thr  <- as.numeric(Fpar$threshold)       # threshold on RDsquare scale

# --- Base R plot ------------------------
eps    <- 1e-12
z_data <- log(pmax(rd2, eps))
z_min  <- min(z_data, na.rm = TRUE)
z_max  <- max(z_data, na.rm = TRUE)
z      <- seq(z_min, z_max, length.out = 4000)

# f_Z(z) = exp(z) * (1/s) * f_F(exp(z)/s)
f_log <- exp(z) * stats::df(exp(z) / s, df1 = df1, df2 = df2) / s

hist(z_data, breaks = 30, freq = FALSE, main = "",
     xlab = "log(RD^2)", ylab = "Density")
lines(z, f_log, lwd = 2)
abline(v = log(thr), col = "red", lty = 2, lwd = 2)


# --- ggplot: histogram of log(RD^2) + mapped F density ------------------------
df_hist  <- data.frame(z = z_data)
df_curve <- data.frame(z = z, dens = f_log)

ggplot(df_hist, aes(z)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "grey80", color = "white") +
  geom_line(data = df_curve, aes(z, dens), linewidth = 1.1) +
  geom_vline(xintercept = log(thr), color = "red", linetype = "dashed") +
  labs(x = "log(RD^2)", y = "Density",
       title = "Histogram of log(RD^2) with mapped F density") +
  theme_minimal()


# --- ggplot by included/excluded -----------------------------------
ind_incld <- RD_obj$ind_incld
lab <- rep("excluded", length(rd2)); lab[ind_incld] <- "included"

df_hist_lab <- data.frame(
  z = z_data,
  observation = factor(lab, levels = c("included","excluded"))
)

ggplot(df_hist_lab, aes(z, fill = observation)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 color = "black", alpha = 0.6, position = "identity") +
  geom_line(data = df_curve, aes(z, 2*dens), color ="#D55E00", linewidth = 1.1, inherit.aes = FALSE) +
  geom_vline(xintercept = log(thr), linetype = "dashed",
             color = "red", linewidth = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
  labs(x = "log(RD^2)", y = "Density",
       title = "log(RD^2) histogram by observation with mapped F density") +
  theme_minimal()


