# visualizations of threshold values using the rrobot outputs directly
# lets load the object from the local folder
abide1_obj <- readRDS("~/Documents/GitHub/RobOutlier-paper/results/ABIDE/ABIDE1.rds")

RD_abide1 <- abide1_obj$RD_obj$RD
RD_sh_abide1 <- abide1_obj$RD_obj_shash$RD

df1 <- data.frame(
  RD = abide1_obj$RD_obj$RD,
  label = ifelse(seq_along(abide1_obj$RD_obj$RD) %in% abide1_obj$RD_obj$ind_incld, "included", "excluded"),
  dataset = "ABIDE1"
)

df1$label <- factor(df1$label, levels = c("included", "excluded"))
