abide1shash_obj <- readRDS("~/Documents/GitHub/RobOutlier-paper/results/ABIDE/ABIDE1_SHASH.rds")
# df_RD <- readRDS("~/Documents/GitHub/RobOutlier-paper/results/df_RD.rds")
# df_HR <- readRDS("~/Documents/GitHub/RobOutlier-paper/results/df_HR.rds")

Zhk <- abide1shash_obj$x_norm

RD_obj_shash <- abide1shash_obj$RD_obj_shash
RDs   <- RD_obj_shash$RD
incl  <- RD_obj_shash$ind_incld
excl  <- RD_obj_shash$ind_excld
Sstar <- RD_obj_shash$S_star
xN    <- abide1shash_obj$x_norm

c(p = RD_obj_shash$p, n = RD_obj_shash$n, h = RD_obj_shash$h)

summary(RDs[incl])           # for inliers, RD² should be ~ O(p), not thousands
summary(RDs[excl])

# Covariance health
ev <- eigen(Sstar, symmetric = TRUE, only.values = TRUE)$values
range(ev); kappa_S <- max(ev)/min(ev); kappa_S

any(!is.finite(xN)); colSums(!is.finite(xN))  # NA/Inf after SHASH?
apply(xN, 2, sd)


Fpar <- abide1shash_obj$thresholds


library(ggplot2)

lab <- rep("excluded", RD_obj$n)
lab[RD_obj$ind_incld] <- "included"

df <- data.frame(
  RD2 = RD_obj$RD,  # squared distances
  observation = factor(lab, levels = c("included","excluded"))
)

# F-density
df1 <- Fpar$df[1]
df2 <- Fpar$df[2]
s   <- 1 / Fpar$scale      # map F → RD²
thr <- Fpar$threshold      # already on RD² scale

# choose an x-range that shows the F curve (otherwise it’s squashed at ~0)
rd2_excl <- RD_obj$RD_excld
x_max <- max(thr * 6, stats::quantile(rd2_excl, 0.95, na.rm = TRUE))

x <- seq(0, x_max, length.out = 100)
f_RD2 <- stats::df(x / s, df1 = df1, df2 = df2) / s
