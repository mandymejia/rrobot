abide1shash_obj <- readRDS("~/Documents/GitHub/RobOutlier-paper/results/ABIDE/ABIDE1_thr03.rds")
# df_RD <- readRDS("~/Documents/GitHub/RobOutlier-paper/results/df_RD.rds")
# df_HR <- readRDS("~/Documents/GitHub/RobOutlier-paper/results/df_HR.rds")

Zhk <- abide1shash_obj$x_norm

RD_obj_shash <- abide1shash_obj$RD_obj_shash
log_RDshash <- log(RD_obj_shash$RD)

df <- data.frame(
  RD2 = RD_obj_shash$RD,
  observation = factor(lab, levels = c("included","excluded"))
)

Fpar <- abide1shash_obj$thresholds
df1 <- Fpar$df[1]
df2 <- Fpar$df[2]
s <- 1 / Fpar$scale
thr <- Fpar$threshold

rd2_excl <- RD_obj_shash$RD_excld
x_max <- max(thr*6, quantile(rd2_excl, 0.95, na.rm=TRUE))
x <- seq(1e-12, x_max, length.out = 400)           # RD² grid (not log)
f_RD2 <- df(x / s, df1, df2) / s

##------------------------------------------------------------------------------
# lets look at the F-density first
df1 <- Fpar$df[1]
df2 <- Fpar$df[2]
s   <- Fpar$scale   # mapping F → RD²

x_max <- Fpar$threshold * 6  # go a few times past the cutoff
x <- seq(0, x_max, length.out = 400)

# F density mapped to RD² scale
f_RD2 <- df(x / s, df1 = df1, df2 = df2) / s

library(ggplot2)
ggplot(data.frame(RD2 = x, density = f_RD2), aes(RD2, density)) +
  geom_line(size = 1) +
  geom_vline(xintercept = Fpar$threshold, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = sprintf("Mapped F density on RD² scale (df1=%d, df2=%.2f)", df1, df2),
       x = expression(RD^2), y = "F density")






library(ggplot2)
ggplot(df, aes(RD2, fill = observation)) +
  geom_histogram(aes(y = ..density..), bins = 60, position = "identity", alpha = 0.8) +
  geom_line(data = data.frame(x=x, y=f_RD2), aes(x, y), inherit.aes=FALSE) +
  geom_vline(xintercept = thr, linetype="dashed", color="red", inherit.aes=FALSE) +
  coord_cartesian(xlim = c(1e-12, x_max)) +
  scale_x_log10() +                                  # <- log axis
  theme_minimal()


library(ggplot2)

df1 <- Fpar$df[1]; df2 <- Fpar$df[2]

s   <- 1 / Fpar$scale
thr <- Fpar$threshold

# use the log10 not e
# Option A: dense grid
x <- seq(0, 200, by = 0.05)
f <- df(x / s, df1 = df1, df2 = df2) / s

library(ggplot2)
ggplot(data.frame(RD2 = log(x), density = f), aes(RD2, density)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = log(thr), linetype = "dashed", color = "red") +
  theme_minimal()

# Option B: stat_function
ggplot(data.frame(RD2 = c(0, 200)), aes(RD2)) +
  stat_function(fun = function(z) df(z / s, df1, df2) / s, n = 5000) +
  geom_vline(xintercept = thr, linetype = "dashed", color = "red") +
  labs(x = expression(RD^2), y = "F density") +
  theme_minimal()


###-----
library(ggplot2)
rd2   <- RD_obj_shash$RD
incl  <- RD_obj_shash$ind_incld
df1   <- Fpar$df[1]; df2 <- Fpar$df[2]
thr   <- as.numeric(Fpar$threshold)

s <- 1 / Fpar$scale
stopifnot(abs(qf(0.99, df1, df2)/Fpar$scale - thr) < 1e-6)  # quick sanity check

# data frame for histogram
lab <- rep("excluded", length(rd2)); lab[incl] <- "included"
df  <- data.frame(RD2 = rd2, observation = factor(lab, c("included","excluded")))

# x-range: show the region where the F curve lives
x_max <- max(thr * 6, quantile(rd2[lab=="excluded"], 0.95, na.rm=TRUE))

ggplot(df, aes(RD2, fill = observation)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60, position = "identity", alpha = 0.75) +
  stat_function(fun = function(z) df(z / s, df1, df2) / s,
                n = 5000, linewidth = 1.1, inherit.aes = FALSE) +
  geom_vline(xintercept = thr, linetype = "dashed", color = "red",
             linewidth = 0.9, inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0, x_max)) +
  scale_fill_manual(values = c("included" = "#009E73", "excluded" = "#56B4E9")) +
  labs(x = expression(RD^2), y = "Density",
       title = "RD\u00B2 histogram (all) with mapped F density (excluded theory)",
       subtitle = sprintf("F_{%d, %.2f} mapped via RD\u00B2 = %.4f × F;  threshold = %.3f",
                          df1, df2, s, thr)) +
  theme_minimal(base_size = 12)


df_ex <- subset(df, observation == "excluded")

ggplot(df_ex, aes(RD2)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "#56B4E9", alpha = 0.75) +
  stat_function(fun = function(z) df(z / s, df1, df2) / s,
                n = 5000, linewidth = 1.1, inherit.aes = FALSE) +
  geom_vline(xintercept = thr, linetype = "dashed", color = "red",
             linewidth = 0.9, inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0, x_max)) +
  labs(x = expression(RD^2), y = "Density",
       title = "RD\u00B2 histogram (excluded only) with mapped F density") +
  theme_minimal(base_size = 12)





# base plot
# inputs
rd2 <- RD_obj_shash$RD
Fpar <- abide1shash_obj$thresholds
df1 <- Fpar$df[1]; df2 <- Fpar$df[2]
s   <- 1 / Fpar$scale               # mapping so RD2= s * F
thr <- as.numeric(Fpar$threshold)

# histogram of log(RD2) with density scale
eps <- 1e-12
z_data <- log(pmax(rd2, eps))
h <- hist(z_data, breaks = 30, freq = FALSE, main = "", xlab = "log(RD^2)")

# build the F density on log-scale: f_Z(z) = exp(z) * (1/s) * f_F(exp(z)/s)
z_min <- min(z_data, na.rm=TRUE); z_max <- max(z_data, na.rm=TRUE)
z <- seq(z_min, z_max, length.out = 2000)
f_log <- exp(z) * df(exp(z) / s, df1 = df1, df2 = df2) / s

# overlay
lines(z, f_log, lwd = 2)
abline(v = log(thr), col = "red", lty = 2, lwd = 2)





library(ggplot2)

eps <- 1e-12
z_data <- log(pmax(RD_obj_shash$RD, eps))
df_hist <- data.frame(z = z_data)

df1 <- Fpar$df[1]; df2 <- Fpar$df[2]
s   <- 1 / Fpar$scale
thr <- as.numeric(Fpar$threshold)

z_min <- min(z_data, na.rm=TRUE)
z_max <- max(z_data, na.rm=TRUE)

curve_df <- data.frame(
  z = seq(z_min, z_max, length.out = 4000)
)

curve_df$dens <- exp(curve_df$z) * df(exp(curve_df$z)/s, df1, df2) / s
# without the histogram

ggplot(df_hist, aes(z)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey80", color = "white") +
  geom_line(data = curve_df, aes(z, dens), linewidth = 1.1) +
  geom_vline(xintercept = log(thr), color = "red", linetype = "dashed") +
  labs(x = "log(RD^2)", y = "density",
       title = "Histogram of log(RD^2) with mapped F density") +
  theme_minimal()


# labels (included / excluded)
ind_incld <- RD_obj_shash$ind_incld
Tn <- length(RD_obj_shash$RD)
lab <- rep("excluded", Tn)
lab[ind_incld] <- "included"

# histogram data on log(RD^2)
eps <- 1e-12
z_data <- log(pmax(RD_obj_shash$RD, eps))
df_hist <- data.frame(
  z = z_data,
  observation = factor(lab, levels = c("included","excluded"))
)

# F density mapped to log(RD^2): f_Z(z) = exp(z) * (1/s) * f_F(exp(z)/s)
df1 <- Fpar$df[1]; df2 <- Fpar$df[2]
s   <- 1 / Fpar$scale
thr <- as.numeric(Fpar$threshold)

z_min <- min(z_data, na.rm = TRUE)
z_max <- max(z_data, na.rm = TRUE)
curve_df <- data.frame(
  z = seq(z_min, z_max, length.out = 4000)
)
curve_df$dens <- exp(curve_df$z) * df(exp(curve_df$z) / s, df1, df2) / s
# install.packages("colorspace")
library(colorspace)
# plot
ggplot(df_hist, aes(z, fill = observation)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, color = "black", alpha = 0.6, position = "identity") +
  geom_line(data = curve_df, aes(x = z, y =2* dens),
            linewidth = 1.1, color = "#D55E00",linetype = 1, inherit.aes = FALSE) +  # <- key
  geom_vline(xintercept = log(thr), linetype = "dashed",
             color = "#D55E00", linewidth = 1, inherit.aes = FALSE) +     # <- key
  scale_fill_manual(values = c("included" = "#009E73", "excluded" = "#D55E00")) +
  labs(x = "log(RD^2)", y = "Density",
       title = "log(RD^2) histogram by observation with mapped F density") +
  theme_minimal()
